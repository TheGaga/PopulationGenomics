# Author:  DINDIN Meryll
# Date:    06 May 2019
# Project: Population Genomics

import os
import time
import tqdm
import joblib
import argparse
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from io import StringIO
from itertools import islice
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import IncrementalPCA