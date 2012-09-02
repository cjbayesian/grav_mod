#!/bin/bash

rsync -rva --delete --exclude-from="server_excludes" ./ darwin:~/CAISN/$1
