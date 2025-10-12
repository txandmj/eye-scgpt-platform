#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

source trapdebug
conda run -n scgpt_cpu python step1_preprocess.py
