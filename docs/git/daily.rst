Consolidating updates
---------------------

When developing your project -- either just adding input files and analysis scripts,
or when doing actual development of new code -- it is important to commit snapshots
to your local working copy of DISPATCH often (essentially after each file edit), but
also important to reorganize those commits before pushing to the on-line repository,
which is typiclally a fork you may be sharing with collaborators.

The dilemma is that pushing all the small incremental commits that you make to the local
working copy generates too much "noise" -- both literally, in the form of emails, and
in the sense that those commits are probably too fine-grained to be useful for others.
On the other hand, refraining to make frequent commits obliviates one of the main advantages
of code maintenance:  the ability to "backtrack", to find out when and what went wrong,
after the fact.

The way out of that dilemma is to reset to the start of the day when done, after saving
the detailed edits as a temporary branch (e.g. named yy-mm-dd for  year-month-day), and
then reorder and clean up at the end of the day (or the end of the week or whenever).
One way of doing that is illusrated here:::

   # -- at start of the day --
   git checkout master                          # get master branch
   git pull                                     # pull in changes from elsewhere
   
   # -- during the day --
   ... edit ...                                 # a small number of edits
   git commit -a -m comment                     # GIT snapshot, with only terse comment
   ... edit ...                                 # a small number of edits
   git commit -a -m comment                     # GIT snapshot, with only terse comment
   ...

   # -- merge in edits from elsewhere --
   git fetch                                    # get updates from elsewhere
   git rebase origin/master                     # your edits become relative to that
   ... may require some edits + commits ...     # in case of collisions
   
   # -- consolidate --
   git branch yy-mm-dd                          # create a backup, temporary branch
   git reset origin/master                      # reset, to get all edits merged together
   git add -p                                   # add _selected_, related edits
   git commit -m "DBG comrehensive comment"     # possibly use 'commit -m' and an editor
   git add -p                                   # add _selected_, related edits
   git commit -m "ORG comrehensive comment"     # possibly use 'commit -m' and an editor
   git add -p                                   # add _selected_, related edits
   git commit -m "DIA comrehensive comment"     # possibly use 'commit -m' and an editor
   git push                                     # push your consolidated commits
   git branch -D yy-mm-dd                       # optionally, delete the temporary branch

One of the advantages with this is that debug statements that were first added and then
removed disappear, reducing the "noise" in the commits that are pushed to the net.

.. toctree::
   :maxdepth: 4

