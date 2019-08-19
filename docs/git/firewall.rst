Behind firewall
----------------

To use GIT with bitbucket on a cluster with a firewall you an use git over SSH:

1. You need to set up bitbucket for SSH (see
   `Atlassians instructions <https://confluence.atlassian.com/bitbucket/set-up-ssh-for-git-728138079.html>`_,
   including the addition of an SSH key to your bitbucket profile.


2. Set-up an SSH tunnel (e.g., see
   `here <https://wickie.hlrs.de/platforms/index.php/Secure_Shell_ssh#Git>`_:

::

  #!bashsession
   
  ssh -R 7777:bitbucket.org:22 username@cluster.domain

or add to your ``$(HOME)/.ssh/config`` on the host / laptop that you use to login to the cluster:

::

  Host clusteralias
      HostName cluster.domain
      User username
      RemoteForward 7777 bitbucket.org:22


3. The first time, when you clone, you need to use an *SSH version* of the bitbucket link, e.g.::

   git clone ssh://git@localhost:7777/aanordlund/dispatch.git

Once cloned, a regular ``git pull`` should work as long as the SSH tunnel is connected.


.. toctree::
   :maxdepth: 4

