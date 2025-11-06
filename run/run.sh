#!/bin/bash
gcc main.c -o main.x -lm -O3

fn="example_subject"
folder="./"

./main.x $fn $folder
