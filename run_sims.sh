#!/bin/bash
conda activate tf
python3 uvn-lhs_10k.py
python3 mvn-lhs_10k.py

python3 uvn-lhs.py
python3 mvn-lhs.py
python3 uvn-rand.py
python3 mvn-rand.py

python3 uvu-lhs-trans.py
python3 uvu-lhs-untr.py
python3 uvu-rand-trans.py
python3 uvu-rand-untr.py