#!/bin/bash

graphlan_annotate.py --annot graphlan-annotation.example.txt graphlan-guide.txt graphlan-annotation.example.xml
graphlan.py graphlan-annotation.example.xml graphlan-annotation.example.png --dpi 300 --size 4.5 --pad 0
