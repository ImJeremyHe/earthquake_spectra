import re


##读取PEER地震动记录文件，返回信息。
def read_record(path):
    ag_list = list()
    dt = 0.0
    record_info = {}
    with open(path, 'r') as f:
        content = f.readlines()
        for line in content[0:4]:
            index = re.search(r'DT=(.*)SE', line)
            if index is not None:
                dt = float(index.group(1).strip())
                break
        for line in content[4:]:
            a = line.split()
            for each in a:
                data = float(each)*9807
                ag_list.append(data)
        #ag_list = [float(x) for x in ag_list]
        ag_list.append(0)
    record_info['ag'] = ag_list
    record_info['dt'] = dt
    return record_info

if __name__ == '__main__':
    info = read_record('1.AT2')
    print(info['ag'])
    print(info['dt'])
