#ssh www-prod mkdir -p /var/www/html/dag.compbio.dundee.ac.uk/share/marek/covid_prot
#scp .htaccess doc/analysis.html www-prod:/var/www/html/dag.compbio.dundee.ac.uk/share/marek/covid_prot


rsync -rvm --include='*/' --include='.htaccess' --include='doc/report.html' --include='tab/*' --exclude='*' . cluster:/cluster/gjb_lab/mgierlinski/public_html/covid_prot