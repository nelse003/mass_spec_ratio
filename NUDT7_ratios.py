import pandas as pd
import os
import re
from io import StringIO

def grouped(iterable, n):
    """
    Chunking iterable into sizes n

    Notes
    ----------
    Python3 uses zip, python 2 would require izip from itertools

    Examples
    ----------
    s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), (s2n,s2n+1,s2n+2,...s3n-1), ...
    """
    return zip(*[iter(iterable)]*n)

def read_grouped_csv(csv):

    """
    Open a csv split by headers preceeed by #

    Ouputs a dictionary of dataframes, prefixed by the header.

    :param csv:
    :return:
    """
    date = os.path.basename(csv).split('_')[0]

    with open(csv) as f:
        data = f.read()

    header_re = re.compile("#*")
    split_data = header_re.split(data)
    split_data.pop(0)

    df_dict = {}
    for header, data in grouped(split_data, 2):
        df_dict[date + " : " + header] = pd.read_csv(StringIO(data))

    return df_dict

if __name__ == "main":
    data_dir = "/hdlocal/enelson/mass_spec_ratio/NUDT7_Data"

    df_dict = {}
    for csv in os.listdir(data_dir):
        csv = os.path.join(data_dir, csv)
        df_dict.update(read_grouped_csv(csv))



