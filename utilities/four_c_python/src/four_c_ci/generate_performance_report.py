# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

import argparse
import dataclasses
import datetime
from typing import Tuple
import os
import json
import dateutil.parser
import plotly.express as px
import jinja2
import pandas as pd


@dataclasses.dataclass
class PerformanceItem:
    min: float
    mean: float
    max: float
    mean_all: float


@dataclasses.dataclass
class TimerResults:
    timestamps: list[datetime.datetime]
    results: list[PerformanceItem]
    commit_hashes: list[str]


@dataclasses.dataclass
class PerformanceResults:
    name: str
    timers: dict[str, TimerResults]


def parse_as_utc(date_str: str) -> datetime.datetime:
    """Parses the given date string as a UTC datetime object."""
    parsed_timestamp = dateutil.parser.parse(date_str)
    if parsed_timestamp.tzinfo is None or parsed_timestamp.utcoffset() is None:
        return parsed_timestamp.replace(tzinfo=datetime.timezone.utc)
    else:
        return parsed_timestamp.astimezone(datetime.timezone.utc)


def read_performance_results(data_file: str) -> list[PerformanceResults]:
    """Reads the performance results from the given JSON file and organizes them into a list of PerformanceResults."""
    with open(data_file, "r") as f:
        all_data = json.load(f)

    performance_tests = {}

    for item in all_data:
        for name, data in item["data"].items():
            if name not in performance_tests:
                performance_tests[name] = PerformanceResults(name, timers={})

            performance_test = performance_tests[name]
            for timer_name in data["Timer names"]:
                if timer_name not in performance_tests[name].timers:
                    performance_tests[name].timers[timer_name] = TimerResults(
                        [], [], []
                    )

                timer = performance_test.timers[timer_name]

                timer.timestamps.append(parse_as_utc(item["date"]))
                timer.commit_hashes.append(item["sha"])
                timer.results.append(
                    PerformanceItem(
                        min=data["Total times"][timer_name]["MinOverProcs"],
                        mean=data["Total times"][timer_name]["MeanOverProcs"],
                        max=data["Total times"][timer_name]["MaxOverProcs"],
                        mean_all=data["Total times"][timer_name]["MeanOverCallCounts"],
                    )
                )

    return list(performance_tests.values())


def generate_plotly_report(
    performance_result: PerformanceResults,
    target_file: str,
    date_min_max: Tuple[datetime.datetime, datetime.datetime],
) -> None:
    """Generates a plotly report for a single performance test and saves it to the target html file."""
    all_data = []
    for timer_name, data in performance_result.timers.items():
        for timestamp, result, sha in zip(
            data.timestamps, data.results, data.commit_hashes
        ):
            all_data.append(
                {
                    "Timestamp": timestamp,
                    "Mean Value": result.mean,
                    "Timer Name": timer_name,
                    "Commit Hash": sha,
                    "Commit URL": f"https://github.com/4C-multiphysics/4C/commit/{sha}",
                }
            )

    df = pd.DataFrame(all_data)
    fig = px.scatter(
        df,
        x="Timestamp",
        y="Mean Value",
        color="Timer Name",
        trendline="ewm",
        trendline_options=dict(halflife=4),
        # trendline_options=dict(window=3),
        title=None,
        hover_data=["Commit Hash", "Commit URL"],
        template="plotly_white",
    )

    fig.update_layout(
        xaxis_title="date of measurement",
        yaxis_title="runtime in seconds",
    )
    fig.update_xaxes(range=[date_min_max[0], date_min_max[1]])

    fig.update_traces(
        marker=dict(
            size=5,
            opacity=0.5,
            symbol="x",
        ),
        line=dict(width=3),
        hovertemplate="<b>%{fullData.name}</b><br />"
        + "Runtime: %{y}s<br />"
        + "Timestamp: %{x}<br />"
        + "<a href='%{customdata[1]}'>%{customdata[0]}</a><extra></extra>",
    )

    visible_timers = ["Input", "Calculation"]
    fig.for_each_trace(
        lambda trace: (
            trace.update(visible="legendonly")
            if trace.name not in visible_timers
            else ()
        )
    )

    fig.write_html(target_file)


def get_date_min_max(
    performance_results: list[PerformanceResults],
) -> tuple[datetime.datetime, datetime.datetime]:
    """Determines the minimum and maximum timestamp across all timers in the performance results."""
    min_date = datetime.datetime.now(datetime.timezone.utc)
    max_date = datetime.datetime.now(datetime.timezone.utc)

    for performance_result in performance_results:
        for timer in performance_result.timers.values():
            min_date = min(min_date, min(timer.timestamps))
            max_date = max(max_date, max(timer.timestamps))

    return min_date, max_date


def main():
    parser = argparse.ArgumentParser(
        prog="generate-performance-report",
        description="Generates a report based on the performance tests.",
    )
    parser.add_argument("test_data")
    parser.add_argument("output_dir")

    args = parser.parse_args()

    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir)

    print(f"Generating report from {args.test_data} to {args.output_dir}...")

    # Read performance results
    performance_results = read_performance_results(args.test_data)

    # determine date time range for x axis (for all tests)
    date_min_max = get_date_min_max(performance_results)

    # generate plotly report for each test
    test_overview = []
    for test in performance_results:

        def sanitize_filename(name: str) -> str:
            return "".join(c if c.isalnum() else "_" for c in name)

        file_name = f"{sanitize_filename(test.name)}.html"
        target_file = os.path.join(args.output_dir, file_name)
        generate_plotly_report(test, target_file, date_min_max)

        # append to overview of reports
        test_overview.append({"title": test.name, "url": file_name})

    # generate overview html file
    env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(
            os.path.join(os.path.dirname(__file__), "resources")
        ),
        autoescape=True,
    )
    template = env.get_template("template.html.j2")
    output_from_parsed_template = template.render(performance_tests=test_overview)
    with open(os.path.join(args.output_dir, "report.html"), "w") as f:
        f.write(output_from_parsed_template)


if __name__ == "__main__":
    main()
