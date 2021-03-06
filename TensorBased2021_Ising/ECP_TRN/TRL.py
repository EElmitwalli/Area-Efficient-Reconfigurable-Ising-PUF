#!/usr/bin/python
import sys, os

import tensorflow as tf
import numpy as np


def ecp_trn(x, y1, rank, n_outputs, ising_size, connection_size):
    weight_initializer = tf.contrib.layers.xavier_initializer(0)
    input_shape = y1.get_shape().as_list()[1:]

    bias = tf.get_variable("bias_{}".format(np.prod(n_outputs)), shape=(1, np.prod(n_outputs)))

    rank1_tnsrs = []

    for i in range(rank):
        rank1_tnsr = []

        for j in range(len(input_shape)):
            rank1_tnsr.append(tf.get_variable("rank1_tnsr_{0}_{1}_{2}".format(i, j, np.prod(n_outputs)),
                                              shape=(input_shape[j]),
                                              initializer=weight_initializer))

        rank1_tnsr.append(tf.get_variable("rank1_tnsr_{0}_output_{1}".format(i, np.prod(n_outputs)),
                                          shape=(n_outputs),
                                          initializer=weight_initializer))

        rank1_tnsrs.append(rank1_tnsr)

    x = tf.reshape(x, [-1, ising_size * ising_size * connection_size])
    cout = tf.zeros([n_outputs], tf.float32)
    for j in range(0, len(rank1_tnsrs)):
        tout = tf.multiply(tf.scalar_mul(1, tf.matmul(x[:, 0:connection_size],
                                                      tf.reshape(rank1_tnsrs[j][0], [-1, 1]))),
                           tf.scalar_mul(1, tf.matmul(
                               x[:, connection_size:connection_size * 2],
                               tf.reshape(rank1_tnsrs[j][1], [-1, 1]))))
        for k in range(2, len(rank1_tnsrs[j]) - 1):
            tout = tf.multiply(tout, tf.scalar_mul(1, tf.matmul(x[:, k*connection_size:(k+1)*connection_size], tf.reshape(rank1_tnsrs[j][k], [-1, 1]))))
        tout = tf.multiply(tout, tf.reshape(rank1_tnsrs[j][k + 1], [-1, 1]))
        cout = tf.add(cout, tout)
    # cout=tf.multiply(tf.scalar_mul(10,wscale),cout)
    return tf.add(cout, bias)

#  mini_x = train_features[j * batch_size:(j + 1) * batch_size, :]
#  x = tf.reshape(x[:,0:connection_size], [-1, connection_size])
#  x = tf.reshape(x[:,k*connection_size:(k+1)*connection_size], [-1, connection_size])
