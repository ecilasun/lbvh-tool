# BVH8 Tool

A tool to generate simple BVH8 structures out of raw triangles


# Summary

This tool will read a list of raw triangle vertices, and generate a BVH8 structure around it. For a BVH8, each node has up to 8 child nodes, and the structure is encoded so that the hardware on E32E can work comfortably with it.

## Building and running on Linux

install 'visual studio code'
sudo apt update
sudo apt upgrade
sudo apt install libx11-dev libsdl2-dev gdb
open the root folder with bhv8-tool in vscode
switch configuration in vscode to 'gbd (Launch)'
use 'configure' from build menu (ctrl+shift+b) (only needed once)
use 'build' from build menu
