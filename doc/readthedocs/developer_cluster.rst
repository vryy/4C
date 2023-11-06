Cloning BACI on a cluster
--------------------------

In general there are several ways to get the Baci source code to a cluster. In the following some possibilities are explained.

**Note:**

- Limitations due to firewall settings may emerge.
- In case you want to commit changes in your repository on a cluster,
  remember to set up Baci as described in the Baci ``README.md`` or in the :ref:`Setup guide <SetupGuidetoBACI>`.
  However, in general, it is recommended to develop code primarily on your workstation.

Method 1: Clone Baci repository from Gitlab
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Clone the Baci repository from Gitlab to a cluster following the procedure as described in ``README.md``.
**Note:** This requires that firewall settings allow access from the cluster to the Gitlab repository.

Method 2: Clone Baci from local workstation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Clone the Baci repository from your workstation to a cluster using ssh::

    ssh <coworker>@<cluster>
    cd <clusterworkspace>
    git clone ssh://<coworker>@<host>/<coworkerworkspace>/baci baci
    exit

**Note:** This requires that firewall settings allow access from the cluster to your workstation.

Method 3: Clone Baci from local workstation via reverse ssh tunnel
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the event that firewall settings forbid direct access from the cluster to your workstation,
open a reverse ssh tunnel from your workstation to the cluster

::

    ssh -R <port>:<workstation>:22 <coworker>@<cluster>

where ``<port>`` is in the range from 1024 through 49151.
> **Note:** Always close the ssh tunnel using the ``exit`` keyword afterwards!

This opens a new shell on the cluster that allows access to your workstation via the ssh tunnel::

    cd <clusterworkspace>
    git clone ssh://<coworker>@localhost:<port>/<coworkerworkspace>/baci baci
    exit


The ssh tunnel access to the workstation is needed every time you want to update anything on the cluster,
e.g., for executing ``git fetch``, ``git pull``, and even ``git push``.
