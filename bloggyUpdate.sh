#!/bin/bash
cd public
(cd ..; $HOME/bin/hugo-0.48/hugo --theme=hugo-geo)
git add --all
echo "Commit message?"
read update
git commit -m "$update"
git push origin gh-pages