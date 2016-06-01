#!/bin/bash
spark-submit --name "$1" --master local[8] --conf spark.eventLog.enabled=true --conf spark.eventLog.dir=/tmp/spark-events digest.py
mv "$1"_* "$2"
