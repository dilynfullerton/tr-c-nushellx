#!/bin/bash

# cougar -> here
from='cougar:~/nushellx/linux/calculations/t0/results/Z2';
to='./results/';

# here -> itheory
# from='./results/Z2';
# to='itheory:~/workspace/tr-c-nushellx/results/';

# rsync -Hrl $from $to
rsync -Hrl \
  --include="*/" \
  --include="*.lpt" \
  --include="*.bat" \
  --include="*.int" \
  --include="*.ans" \
  --exclude="*" \
  $from $to;
