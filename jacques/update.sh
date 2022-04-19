#!/bin/bash

# Update JACQUES and DYNAMON software from GitHub

# ask for confirmation
read -p "Update JACQUES and DYNAMON software? (y/n): " response

if [ "${response,,}" == "y" -o "${response,,}" == "yes" ]; then
  # clean all compiled files of DYNAMON and fDynamo
  printf "$(tput setaf 5)Cleaning $(tput setaf 3; tput bold)DYNAMON$(tput sgr0) \n"
  make clean_all -C $JACQUES/dynamon/src

  # update JACQUES
  printf "$(tput setaf 5)Updating $(tput setaf 3; tput bold)JACQUES$(tput sgr0) \n"
  git -C $JACQUES fetch --all
  git -C $JACQUES reset --hard origin/master

  # update DYNAMON
  printf "$(tput setaf 5)Updating $(tput setaf 3; tput bold)DYNAMON$(tput sgr0) \n"
  git -C $JACQUES submodule update --init --recursive

  # compile DYNAMON
  printf "$(tput setaf 5)Compiling $(tput setaf 3; tput bold)DYNAMON$(tput sgr0) \n"
  make -C $JACQUES/dynamon/src
else
  printf "Update cancelled. \n"
  exit 1
fi
