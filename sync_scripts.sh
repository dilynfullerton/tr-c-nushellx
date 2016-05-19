#!/bin/bash
# Put the python scripts, as they are in this directory, onto cougar.
from='./*.py';
to='cougar:~/nushellx/linux/calculations/t0/';
rsync -r ${from} ${to};
