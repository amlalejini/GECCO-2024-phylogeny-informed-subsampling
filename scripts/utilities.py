import errno, os, csv, copy

def mkdir_p(path):
    """
    This is functionally equivalent to the mkdir -p [fname] bash command
    """
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

def extract_params_cmd_log(path, exec_name="diagnostics-suite"):
    """
    Extract Avida parameters from log of command used to run Avida.
    The log should contain only the text used to run Avida.

    e.g. a cmd log containing,
        ./diagnostics-suite -DIAGNOSTIC struct-exploitation -POP_SIZE 500 -SEED 1055 -SELECTION lexicase -MAX_GENS 50000
    """
    content = None
    with open(path, "r") as fp:
        content = fp.read().strip()
    content = content.replace(f"./{exec_name}", "").strip()
    cmd_args = content.split(" ")
    cfg = {}
    for i in range(0, len(cmd_args), 2):
        assert i+1 < len(cmd_args)
        key = cmd_args[i].strip("-")
        cfg[key] = cmd_args[i+1]
    return cfg

def read_csv(file_path):
    content = None
    with open(file_path, "r") as fp:
        content = fp.read().strip().split("\n")
    header = content[0].split(",")
    content = content[1:]
    lines = [{header[i]: l[i] for i in range(len(header))} for l in csv.reader(content, quotechar='"', delimiter=',', quoting=csv.QUOTE_ALL, skipinitialspace=True)]
    return lines

# NOTE - this function assumes that data will be ordered!
def filter_ordered_data(data, units, resolution):
    if (units == "interval"):
        return filter_ordered_data_interval(data, resolution)
    elif (units == "total"):
        return filter_ordered_data_total(data, resolution)
    elif (units == "evals"):
        return filter_ordered_data_keyval(data, "evals", resolution)
    elif (units == "gens"):
        return filter_ordered_data_keyval(data, "gens", resolution)
    else:
        return data

def filter_ordered_data_keyval(data, key, interval):
    prev_val = 0
    ret_data = []
    for i in range(len(data)):
        cur_val = int(data[i][key])
        if (i==0) or cur_val >= (prev_val+interval) or (i==(len(data)-1)):
            ret_data.append(copy.deepcopy(data[i]))
            prev_val = cur_val
    return ret_data

def filter_ordered_data_interval(data, interval):
    return [data[i] for i in range(len(data)) if (i==0) or (not (i % interval)) or (i==len(data)-1)]

def filter_ordered_data_total(data, total):
    sample = [int(x*(len(data)-1)/(total-1)) for x in range(total)]
    return [data[i] for i in sample]
