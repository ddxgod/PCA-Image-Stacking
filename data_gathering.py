import zipline as zp
import quandl as qd
import pandas as pd
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import pytz
import os.path


start_datestart_d  = '2015-01-01'
end_date = '2018-01-01'

companies = ['AAPL', 'MSFT', 'AMZN', 'INTC', 'TSLA']
# Grab data for select companies
qd.ApiConfig.api_key = 'Zx2BTzzcz264ssx748NF'
data = qd.get_table('WIKI/PRICES', ticker = companies, date = {'gte': start_date, 'lte': end_date}, paginate=True)

data.head()