#!/bin/bash
chown shiny:shiny -R /srv/shiny-server/digest
sudo -u shiny bash << EOF
shiny-server&
EOF

