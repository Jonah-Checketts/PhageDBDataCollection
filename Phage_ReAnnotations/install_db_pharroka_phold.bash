#!/bin/bash
#Do this second (dont have to rerun it after it has been ran once)

echo "Downloading pharokka database. This will take a few minutes. Please be patient :)"
time install_databases.py -o pharokka_db

echo "Downloading phold database. This will take a few minutes. Please be patient :)"
time phold install -d phold_db