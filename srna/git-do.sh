#! /usr/bin/env bash
set -eu -o pipefail
_wd=$(pwd)
_path=$(dirname $0 | xargs -i readlink -f {})

#### git clone
git clone https://github.com/d2jvkpn/srna
cd srna

#### setup
git remote set-url origin git@github.com:d2jvkpn/srna.git
git remote set-url --add origin git@gitlab.com:d2jvkpn/srna.git
git config user.name d2jvkpn
git config user.email chenbin2018@protonmail.com

git branch -C dev
git checkout dev

git commit -am "branch dev"

git push --set-upstream origin dev
git branch -a

#### commit
git add -A
git commit -m "init"
git push
