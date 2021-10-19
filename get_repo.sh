#!/bin/bash
## Script that clones this repository
cd /home/jovyan/work

# git reflog requires a name and email if user is not in passwd
# even if you're only cloning
export GIT_COMMITTER_NAME=anonymous
export GIT_COMMITTER_EMAIL=anon@localhost

git clone git@github.com:csc-training/r-spatial-course.git
