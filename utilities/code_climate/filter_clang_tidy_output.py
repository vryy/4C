import re
import sys

begin_warning_regex = re.compile(r"^.* warning: .*\[.*\]$")
begin_note_regex = re.compile(r"^\/\S+\:\d+\:\d+\: note: .*$")
note_regex = re.compile(r"^note: .*$")
message_block_regex = re.compile(r"^ *\d* \|.*$")
color_output_regex = re.compile(r"\033\[[0-9;]+m")

is_warning_or_note_active = False
for line in sys.stdin:
    line_removed_coloroutput = color_output_regex.sub("", line)

    # skip all warning message blocks
    if is_warning_or_note_active and message_block_regex.match(
        line_removed_coloroutput
    ):
        continue

    is_warning_or_note_active = False

    # skip warnings
    if begin_warning_regex.match(line_removed_coloroutput) or begin_note_regex.match(
        line_removed_coloroutput
    ):
        is_warning_or_note_active = True
        continue

    # ignore pure notes
    if note_regex.match(line_removed_coloroutput):
        continue

    # write all other lines
    sys.stdout.write(line)
