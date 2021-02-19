#!/bin/bash

# Description: Update JACQUES software from GitHub
# Last update: 28-12-2020

# ask for confirmation
read -p "Update JACQUES and DYNAMON software? (y/n): " response

if [ "${response,,}" == "y" ] || [ "${response,,}" == "yes" ]; then
  workdir=$PWD
  cd $JACQUES

  # clean all compiled files of DYNAMON and fDynamo
  printf "$(tput setaf 5)Cleaning $(tput setaf 3; tput bold)DYNAMON$(tput sgr0) \n"
  make clean_all -C dynamon/src

  # update JACQUES
  printf "$(tput setaf 5)Updating $(tput setaf 3; tput bold)JACQUES$(tput sgr0) \n"
  git fetch --all
  git reset --hard origin/master

  # update DYNAMON
  printf "$(tput setaf 5)Updating $(tput setaf 3; tput bold)DYNAMON$(tput sgr0) \n"
  git submodule update --init --recursive

  # compile DYNAMON
  printf "$(tput setaf 5)Compiling $(tput setaf 3; tput bold)DYNAMON$(tput sgr0) \n"
  make -C dynamon/src

  cd $workdir
else
  printf "Update cancelled. \n"
  exit
fi
