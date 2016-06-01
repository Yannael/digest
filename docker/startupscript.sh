#!/bin/bash
#/etc/init.d/openvpn restart
chown shiny:shiny -R /srv/shiny-server/digest
sudo -u shiny bash << EOF
shiny-server&
ssh -nNT -p 23 -L 21050:127.0.0.1:21050 digest@hadoop.mlg.ulb.ac.be&
EOF

