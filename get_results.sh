#!/bin/bash
from='cougar:~/nushellx/linux/calculations/t0/results/Z2';
to='./results/';
# rsync -Hrl $from $to
rsync -Hrl \
  --include="*/" \
  --include="*.lpt" \
  --include="*.bat" \
  --include="*.int" \
  --include="*.ans" \
  --exclude="*" \
  $from $to;
