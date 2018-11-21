import tensorflow as tf
import numpy as np
import os, sys
import time
import math

# global parameters
resolution = 100

class Param:
    def __init__(self, label, struc, start, end):
        self.label = label
        self.struc = struc
        self.start = start
        self.end   = end

valid_param = Param("../../fenics/run_1_valid", "../../utility/run_1_valid/output", 0, 99)
train_param = Param("../../fenics/run_1_valid", "../../utility/run_1_valid/output", 100, 499)

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

def get_struc(path, index):
    '''
    input:
    path + index[i] defines a stucture file
    output:
    (len(index), dim_z, dim_y, dim_x, 1)
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
        struc.append(data)
    return np.array(struc)[:, :, :, :, np.newaxis]

def get_tf_input(labelpath, structpath, start, end):
    '''
    input:
    path labelpath containing result.txt and structpath containing 3D_n.dat
    start and end defines index range
    output:
    label: (end-start+1,)
    struc: (end-start+1, dim, dim, dim, 1)
    '''
    index, label = get_label(labelpath + "/result.txt", start, end)
    struc = get_struc(structpath, index)
    return struc, label

def get_all_data():
    '''
    get training data and validating data
    '''
    # validating data
    log_file = open("log.txt", "a")
    time_start=time.time()
    valid_struc, valid_label = get_tf_input(valid_param.label, valid_param.struc, valid_param.start, valid_param.end)
    time_end=time.time()
    print("\nX_valid: ", valid_struc.shape, file = log_file)
    print("y_valid: ", valid_label.shape, file = log_file)
    print(time_end-time_start, file = log_file)
    log_file.close()
    # training data
    log_file = open("log.txt", "a")
    time_start=time.time()
    train_struc, train_label = get_tf_input(train_param.label, train_param.struc, train_param.start, train_param.end)
    time_end=time.time()
    print("\nX_train: ", train_struc.shape, file = log_file)
    print("y_train: ", train_label.shape, file = log_file)
    print(time_end-time_start, file = log_file)
    log_file = open("log.txt", "a")
    # normalize to mean 0, std 1
    mean, std = np.mean(train_label), np.std(train_label)
    train_label = (train_label - mean) / std
    valid_label = (valid_label - mean) / std
    # save results
    log_file = open("log.txt", "a")
    print("mean = {0:.3g}, std = {1:.3g}".format(mean, std), file = log_file)
    log_file.close()
    np.savetxt("y_valid.txt", valid_label, fmt = '%.4f')
    np.savetxt("y_train.txt", train_label, fmt = '%.4f')
    return valid_struc, valid_label, train_struc, train_label

def complex_model(X, is_training):
    # define weight
    Wconv1 = tf.get_variable("Wconv1", shape=[3, 3, 3, 1, 16])
    bconv1 = tf.get_variable("bconv1", shape=[16])
    Wconv2 = tf.get_variable("Wconv2", shape=[3, 3, 3, 16, 32])
    bconv2 = tf.get_variable("bconv2", shape=[32])
    Wconv3 = tf.get_variable("Wconv3", shape=[3, 3, 3, 32, 64])
    bconv3 = tf.get_variable("bconv3", shape=[64])
    # fc1
    Wfc1   = tf.get_variable("Wfc1", shape=[4*4*4*64,64])
    # fc2
    Wfc2   = tf.get_variable("Wfc2", shape=[64,1])
    
    #define graph
    conv1 = tf.nn.conv3d(X,Wconv1,[1,2,2,2,1],padding="SAME") + bconv1
    conv1_norm = tf.layers.batch_normalization(conv1, training=is_training)
    ## 100->50
    relu1 = tf.nn.relu(conv1_norm)
    pool1 = tf.nn.max_pool3d(relu1,[1,2,2,2,1],padding="SAME",strides=[1,2,2,2,1])
    ## 50->25
    
    conv2 = tf.nn.conv3d(pool1,Wconv2,[1,2,2,2,1],padding="SAME") + bconv2
    conv2_norm = tf.layers.batch_normalization(conv2, training=is_training)
    ## 25->13
    relu2 = tf.nn.relu(conv2_norm)
    pool2 = tf.nn.max_pool3d(relu2,[1,2,2,2,1],padding="SAME",strides=[1,2,2,2,1])
    ## 13->7
    
    conv3 = tf.nn.conv3d(pool2,Wconv3,[1,2,2,2,1],padding="SAME") + bconv3
    conv3_norm = tf.layers.batch_normalization(conv3, training=is_training)
    ## 7->4
    relu3 = tf.nn.relu(conv3_norm)

    flat  = tf.reshape(relu3,[-1,4*4*4*64])
    fc1   = tf.matmul(flat,Wfc1)
    reluf = tf.nn.relu(fc1)
    fc2   = tf.matmul(reluf,Wfc2)
    
    return fc2[:, 0]

def run_valid(sess, variables, feed_dict):
    loss, y_pred = sess.run(variables, feed_dict=feed_dict)
    np.savetxt("y_valid_pred.txt", y_pred, fmt = '%.4f')
    return loss

def run_train(sess, X, y, variables_train, feed_dict_train, variables_valid, feed_dict_valid, saver, epochs=1, batch_size=100):
    log_file = open("log.txt", "a")
    los_file = open("los.txt", "a")
    print("\nTraining", file = log_file)
    print("epoch batch     time     loss", file = log_file)
    # shuffle indicies
    Xd, yd = feed_dict_train[X], feed_dict_train[y]
    train_indicies = np.arange(Xd.shape[0])
    np.random.shuffle(train_indicies)
    np.savetxt("y_train.txt", yd[train_indicies], fmt = '%.4f')
    # iterate over epoches
    time_start=time.time()
    for e in range(epochs):
        # keep track of losses and accuracy
        loss_trai = 0
        y_pred_trai = []
        for i in range(int(math.ceil(Xd.shape[0]/batch_size))):
            # generate indicies for the batch
            start_idx = (i*batch_size)%Xd.shape[0]
            idx = train_indicies[start_idx:start_idx+batch_size]
            # create a feed dictionary for this batch
            feed_dict_train[X], feed_dict_train[y] = Xd[idx,:], yd[idx]
            # get batch size
            actual_batch_size = yd[idx].shape[0]
            # have tensorflow compute losses and predictions and do gradient decent
            loss, y_pred, _ = sess.run(variables_train,feed_dict=feed_dict_train)
            loss_trai += loss*loss*actual_batch_size
            y_pred_trai.extend(list(y_pred))
            time_end=time.time()
            print("{0:5d} {1:5d} {2:8.1f} {3:7.3f}".format(e, i, time_end-time_start, loss), file = log_file)
            log_file.close()
            log_file = open("log.txt", "a")
        loss_trai = np.sqrt(loss_trai/Xd.shape[0])
        # testing results
        loss_test, y_pred_test = sess.run(variables_valid, feed_dict=feed_dict_valid)
        # save training and testing results
        print("Epoch {0:5d} training loss = {1:7.3f} testing loss = {2:7.3f}".format(e + 1, loss_trai, loss_test), file = los_file)
        np.savetxt("y_train_pred.txt", y_pred_trai, fmt = '%.4f')
        np.savetxt("y_valid_pred.txt", y_pred_test, fmt = '%.4f')
        # save the model
        if (e % 10 == 0):
            saver.save(sess, "Model/model.ckpt"+str(e))
            log_file.close()
            los_file.close()
            log_file = open("log.txt", "a")
            los_file = open("los.txt", "a")
    log_file.close()
    los_file.close()

def main():
    # testing data and training data
    valid_struc, valid_label, train_struc, train_label = get_all_data()
    # placeholders
    X = tf.placeholder(tf.float32, [None, resolution, resolution, resolution, 1])
    y = tf.placeholder(tf.float32, [None])
    is_training = tf.placeholder(tf.bool)
    # graph
    y_pred = complex_model(X, is_training)
    # loss
    loss = tf.sqrt(tf.losses.mean_squared_error(y_pred, y))
    # optimizer
    optimizer = tf.train.RMSPropOptimizer(1e-3)
    update_ops = tf.get_collection(tf.GraphKeys.UPDATE_OPS)
    with tf.control_dependencies(update_ops):
        training = optimizer.minimize(loss)
    # saver
    saver = tf.train.Saver()
    # run tensorflow
    variables_valid = [loss, y_pred]
    feed_dict_valid = {X: valid_struc, y: valid_label, is_training: True}
    variables_train = [loss, y_pred, training]
    feed_dict_train = {X: train_struc, y: train_label, is_training: True}
    with tf.Session() as sess:
        # restore weights or initialize weights
        # saver.restore(sess, "./Model/model.ckpt" + str(50))
        sess.run(tf.global_variables_initializer())
        # check the initial loss
        log_file = open("log.txt", "a")
        print("\nvalid loss:", run_valid(sess, variables_valid, feed_dict_valid), file = log_file)
        log_file.close()
        # train the model
        run_train(sess, X, y, variables_train, feed_dict_train, variables_valid, feed_dict_valid, saver, 100, 20)

if __name__ == '__main__':
    sys.exit(main())
