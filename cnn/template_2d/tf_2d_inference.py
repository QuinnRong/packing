import tensorflow as tf
import re, os, sys
import time
import numpy as np

# global parameters
resolution = 100
model_file = "./Model/model.ckpt0"

class Param:
    def __init__(self, label, struc, start, end):
        self.label = label
        self.struc = struc
        self.start = start
        self.end   = end

test_param = Param("../../fenics/run_1_valid", "../../utility/run_1_valid/output", 0, 199)

def get_label(filename, start, end):
    '''
    input: txt file, range from start to end
    idx1 val1
    idx2 val2
    ...
    output: an array of idx(string), an array of val(float)
    '''
    index = []
    label = []
    with open(filename) as f:
        lines = f.readlines()
        if len(lines) <= end:
            sys.exit("index out of range!")
        for i in range(start, end + 1):
            index.append(lines[i].split()[0])
            label.append(float(lines[i].split()[1]))
    return np.array(index), np.array(label)

def get_cross_section(dire, dist, struct):
    '''
    input:
    direction: x, y or z
    distance: coordinate in the direction
    struct: array of shape (dim_z, dim_y, dim_x)
    output:
    2 dimensional cross section
    '''
    if dire == "x":
        return struct[:,:,dist]
    elif dire == "y":
        return struct[:,dist,:]
    elif dire == "z":
        return struct[dist,:,:]
    else:
        sys.exit("direction not correct!")

def get_struc(path, index, dire, dist):
    '''
    input:
    path + index[i] defines a stucture file
    dire[i] + dist[i] represent a cross section
    output:
    (len(index), dim_h, dim_w, len(dire))
    '''
    struc = []
    for id in index:
        filename = "3D_" + id + ".dat"
        data = np.loadtxt(os.path.join(path, filename), dtype=int)
        dim_x = data.shape[1]
        if dim_x != resolution:
            sys.exit("resolution not match!")
        if data.shape[0] != dim_x**2:
            sys.exit("input is not a cubic!")
        data = data.reshape((resolution, resolution, resolution))
        image_data = []
        for i in range(len(dire)):
            image_data.append(get_cross_section(dire[i], dist[i], data))
        struc.append(np.array(image_data).transpose(1, 2, 0))
    return np.array(struc)

def get_tf_input(labelpath, structpath, start, end, dire, dist):
    '''
    input:
    path labelpath containing result.txt and structpath containing 3D_n.dat
    start and end defines index range
    output:
    label: (end-start+1,)
    struc: (end-start+1, dim, dim, len(dire))
    '''
    index, label = get_label(labelpath + "/result.txt", start, end)
    struc = get_struc(structpath, index, dire, dist)
    return struc, label

def get_mean_std(filename, nline):
    '''
    input:
    nth line in the file is, for example,
    mean = 2.5549, std = 0.871817
    output:
    float(2.5549), float(0.871817)
    '''
    with open(filename) as f:
        for i in range(nline - 1):
            f.readline()
        line = f.readline()
        linelist = re.split("[=,]", line)
    return float(linelist[1]), float(linelist[3])

def get_all_data():
    '''
    get testing data
    '''
    dire = ["x", "y", "x", "y"]
    dist = [0, 0, 50, 50]
    # testing data
    time_start=time.time()
    test_struc, test_label = get_tf_input(test_param.label, test_param.struc, test_param.start, test_param.end, dire, dist)
    time_end=time.time()
    print("\nX_test: ", test_struc.shape)
    print("y_test: ", test_label.shape)
    print(time_end-time_start)
    mean, std = get_mean_std("log.txt", 9)
    test_label = (test_label - mean) / std 
    np.savetxt("y_test.txt", test_label, fmt = '%.4f')
    return test_struc, test_label

def run_test(sess, variables, feed_dict):
    loss, y_pred = sess.run(variables, feed_dict=feed_dict)
    np.savetxt("y_test_pred.txt", y_pred, fmt = '%.4f')
    print(loss)
    return loss

def main():
    # testing data
    test_struc, test_label = get_all_data()
    # run tensorflow
    with tf.Session() as sess:
        saver = tf.train.import_meta_graph(model_file + ".meta")
        saver.restore(sess, model_file)

        graph = tf.get_default_graph()
        loss = graph.get_operation_by_name("loss").outputs[0]                   # loss:0
        y_pred = graph.get_operation_by_name("y_pred").outputs[0]               # y_pred:0
        X = graph.get_operation_by_name("X").outputs[0]                         # X:0
        y = graph.get_operation_by_name("y").outputs[0]                         # y:0
        is_training = graph.get_operation_by_name("is_training").outputs[0]     # is_training:0
        keep_prob = graph.get_operation_by_name("keep_prob").outputs[0]         # keep_prob:0

        variables_test = [loss, y_pred[:,0]]
        feed_dict_test = {X: test_struc, y: test_label, is_training: True, keep_prob: 1}

        run_test(sess, variables_test, feed_dict_test)

if __name__ == '__main__':
    sys.exit(main())
