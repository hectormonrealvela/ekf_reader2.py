import struct
import numpy as np
import sys, time
import binascii
import octomap
import threading
import time
from pcapfile import savefile
from dynamic import Dynamic
from collections import deque
import matplotlib.pyplot as plt
from transforms3d.euler import euler2mat, mat2euler
import transforms3d.derivations.eulerangles
import pyproj
import pynmea2
import velodyne
import json
from scipy.signal import medfilt, medfilt2d, wiener, spline_filter, cubic
import utm
import datetime
from mpldatacursor import datacursor
import matplotlib


class Thread_EKF(threading.Thread):
    def __init__(self, threadName, ekf_reader):
        threading.Thread.__init__(self)
        self.ekf_reader = ekf_reader
        self.threadName = threadName
        print "Thread created: " + self.threadName

    def run(self):
        print "Starting treading: " + self.threadName
        self.ekf_reader.read()

    def stop(self):
        print "Stoping treading: " + self.threadName
        sys.exit()


class ekf_reader:
    def __init__(self, ekf_file, pcap_file=None, verbose=False, show_data_flag=True):

        self.ACCELEROMETER_SENSITIVITY = 16384.0  # G
        self.GYROSCOPE_SENSITIVITY = 131.0  # dps

        self.file_p = open('file_vuleta_1.txt', 'w')
        self.file_descriptor = open('ekf_file.txt', "r+")
        line = self.file_descriptor.readline()

        if verbose == True:
            self.line_list = line.split()
            for index in range(len(self.line_list)):
                print str(index) + ' ' + str(self.line_list[index])

        self.pcap_file_name = pcap_file
        if self.pcap_file_name != None:
            pointer_pcap_file = open(self.pcap_file_name, 'rb')
            self.pcap_file = savefile.load_savefile(pointer_pcap_file, verbose=False)
            self.index_pkt = 188701
            self.velodyne = velodyne.HDL_32E()

        # sys.exit()
        self.total_duration_gps = 0
        self.init_time_stamp = 0
        self.gps_init_time = 0
        self.gps_time_ms = 0
        self.gps_time_ms_prev = 0
        self.gps_time_ms_prev_aux = 0
        self.gps_time_corrected_ms = []
        self.end_time = 0
        self.line_list = []
        self.follow_line_list = []
        self.verbose = verbose
        self.time = 0
        self.gps_time = 0
        self.data = {}
        self.data_queue = deque([])
        self.is_finished = False
        self.init_gyr_x = 0
        self.gyr_x = 0
        self.gyr_x_list = []
        self.init_gyr_y = 0
        self.gyr_y = 0
        self.gyr_y_list = []
        self.init_gyr_z = 0
        self.gyr_z = 0
        self.gyr_z_list = []
        self.acc_x = 0
        self.acc_x_list = []
        self.acc_y = 0
        self.acc_y_list = []
        self.acc_z = 0
        self.acc_z_list = []
        self.yaw = 0
        self.yaw_offset = 0
        self.yaw_corrected = 0
        self.yaw_prev = 0
        self.yaw_prev_2 = 0
        self.yaw_prev_3 = 0
        self.yaw_prev_4 = 0
        self.yaw_prev_5 = 0
        self.yaw_list = []
        self.yaw_list_corrected = []
        self.yaw_evaluated = 0
        self.yaw_evaluated_prev = 0
        self.yaw_evaluated_list = []
        self.pitch = 0
        self.pitch_list = []
        self.roll = 0
        self.roll_list = []
        self.y_calidad_list = []
        self.x_calidad_list = []
        self.y_list = []
        self.x_list = []
        self.time_stamp_ms = []
        self.init_yaw = 0
        self.init_pitch = 0
        self.init_roll = 0
        self.init_x = 0
        self.init_y = 0
        self.init_x_ekf = 0
        self.init_y_ekf = 0
        self.init_azim = 0
        self.x = 0
        self.y = 0
        self.x_ekf = 0
        self.y_ekf = 0
        self.p1 = pyproj.Proj(init='epsg:32630')
        self.p1 = pyproj.Proj(proj='utm', zone=30)
        self.sum_azim = 0
        self.show_data_flag = show_data_flag
        self.calidad_gps = 0
        self.calidad_gps_list = []

        self.y_ekf_list = []
        self.x_ekf_list = []
        self.time_diff = []

        self.prev_time = 0

        self.x_time_aux = np.array([], dtype=np.float32)
        self.time_inter = np.array([], dtype=np.float32)
        self.x_inter = np.array([], dtype=np.float32)
        self.y_inter = np.array([], dtype=np.float32)
        self.yaw_inter = np.array([], dtype=np.float32)
        self.rol_inter = np.array([], dtype=np.float32)
        self.pit_inter = np.array([], dtype=np.float32)

        self.count = 0
        self.offset = 0
        self.correction = 0
        self.correction_prev = 0

        self.follow_line = None

    def read(self):
        while self.is_finished == False:

            stop_limit_time = 57848000
            if self.pcap_file_name != None:
                pkt = self.pcap_file.packets[self.index_pkt]
                binary_pkt = pkt.raw()
                self.index_pkt += 1

                if len(binary_pkt) == 554:  # GPS data
                    self.velodyne.evaluation(binary_pkt)
                    print self.velodyne.get_timeStamp()

            line = self.file_descriptor.readline()
            last_pos = self.file_descriptor.tell()

            for index in range(1):
                line_new = self.file_descriptor.readline()
                line_new = self.file_descriptor.readline()

            self.file_descriptor.seek(last_pos)

            self.count += 1

            if len(line) != 0 and self.gps_time_ms < stop_limit_time:
                self.line_list = line.split()
                self.follow_line_list = line_new.split()
                self.time = int(self.line_list[0])
                self.gps_time = int(self.line_list[1])
                output = utm.from_latlon(pynmea2.dm_to_sd(self.line_list[2]), pynmea2.dm_to_sd(self.line_list[3]))
                self.x = output[0]
                self.y = output[1]
                self.x = -self.x
                self.x_ekf, self.y_ekf = float(self.line_list[18]), float(self.line_list[17])

                if self.init_time_stamp == 0:
                    EKF_NAntena = 4484750;
                    EKF_EAntena = 470578;

                    self.init_time_stamp = int(self.line_list[0])  # Sampling time 10ms
                    self.gps_init_time = int(self.line_list[1])  # Sampling time 50ms
                    self.init_yaw = float(self.line_list[19])
                    self.init_pitch = float(self.line_list[70])
                    self.init_roll = float(self.line_list[71])
                    self.init_gyr_x = float(self.line_list[67])
                    self.init_gyr_y = float(self.line_list[68])
                    self.init_gyr_z = float(self.line_list[69])
                    print self.x
                    print self.y
                    self.init_x = self.x
                    self.init_y = self.y
                    self.init_x_ekf = EKF_EAntena
                    self.init_y_ekf = EKF_NAntena

                # Time
                time_aux = self.time
                self.gps_time_ms = int(self.line_list[1])

                if self.gps_time_ms != self.gps_time_ms_prev_aux:
                    self.gps_time_corrected_ms.append(self.gps_time_ms)
                    self.gps_time_ms_prev_aux = self.gps_time_ms

                    #  if self.gps_time_ms != self.gps_time_ms_prev_aux:
                    #      self.init_time_stamp = int(self.line_list[0])
                    #      self.gps_time_ms_prev_aux = self.gps_time_ms

                    #      #if self.time_stamp_ms[-1]%50 < 45.0:
                    #      #    self.gps_time_ms_prev_aux = self.gps_time_ms

                    #  self.time_stamp_ms.append((int(self.line_list[0])-self.init_time_stamp)/1000.0)
                    #  utm_hour = time.gmtime(self.gps_time_ms/1000)

                    #  #print str(self.time_stamp_ms[-1]%50) + ' ' + str(self.gps_time_ms_prev_aux) + ' ' + str(self.gps_time_ms_prev_aux + self.time_stamp_ms[-1]%50)
                    #  print str(str(self.gps_time_ms_prev_aux + self.time_stamp_ms[-1]))
                    #  #self.gps_time_corrected_ms.append(self.gps_time_ms_prev_aux + self.time_stamp_ms[-1]%50)
                    #  self.gps_time_corrected_ms.append(self.gps_time_ms_prev_aux + self.time_stamp_ms[-1])

                    #  if self.gps_time_ms != self.gps_time_ms_prev:
                    #     self.offset = ((time_aux%50000)/1000.0) - ((time_aux%50000)/1000.0)%10
                    #     if self.correction < 0:
                    #         self.correction += 50
                    #     self.correction = ((time_aux%50000)/1000.0) - self.offset
                    #     self.gps_time_corrected_ms.append(int(self.line_list[1]) + self.correction )
                    #     self.correction_prev = 0

                    #  else:
                    #     self.correction = ((time_aux%50000)/1000.0) - self.offset
                    #     if self.correction < 0:
                    #         self.correction += 50

                    #     if self.correction < self.correction_prev:
                    #         self.correction += 50

                    #     self.gps_time_corrected_ms.append(int(self.line_list[1]) + self.correction )

                    if self.verbose == True and self.init_time_stamp != 0:
                        print 'Offset: ' + str(self.offset)
                        print 'Correction: ' + str(self.correction)

                    self.gps_time_ms_prev = self.gps_time_ms
                    self.correction_prev = self.correction

                    # Yaw
                    self.yaw = float(self.line_list[19])
                    if abs(self.yaw_prev - self.yaw) > np.pi / 2:
                        if self.yaw_prev > np.pi:
                            if self.yaw_offset == 0:
                                self.yaw_offset = (int(self.yaw_prev / np.pi) + 1) * np.pi
                            else:
                                self.yaw_offset += 2 * np.pi

                    self.yaw_corrected = self.yaw + self.yaw_offset

                    self.yaw_list.append(self.yaw)
                    self.yaw_list_corrected.append(self.yaw_corrected)

                    self.yaw_prev_5 = self.yaw_prev_4
                    self.yaw_prev_4 = self.yaw_prev_3
                    self.yaw_prev_3 = self.yaw_prev_2
                    self.yaw_prev_2 = self.yaw_prev
                    self.yaw_prev = self.yaw

                    rotation_matix = euler2mat(0, 0, 0, 'rxyz')
                    rotation_matix_ekf = euler2mat(0, 0, 0, 'rxyz')

                    result = np.dot([self.x - self.init_x, self.y - self.init_y, 0], rotation_matix)
                    result_ekf = np.dot([self.x_ekf - self.init_x_ekf, self.y_ekf - self.init_y_ekf, 0],
                                        rotation_matix_ekf)

                    # Calidad GPS
                    self.calidad_gps = float(self.line_list[6])
                    self.calidad_gps_list.append(self.calidad_gps)

                    if self.calidad_gps != 4:
                        self.x_calidad_list.append(None)
                        self.y_calidad_list.append(None)
                    else:
                        self.x_calidad_list.append(result_ekf[0])
                        self.y_calidad_list.append(result_ekf[1])

                    self.x_list.append(result[0])
                    self.y_list.append(result[1])


                    self.x_ekf_list.append(result_ekf[0])
                    self.y_ekf_list.append(result_ekf[1])

                    # Yaw evaluation
                    if len(self.x_ekf_list) > 2 and (
                                    self.y_ekf_list[-1] != self.y_ekf_list[-2] or self.x_ekf_list[-1] !=
                                self.x_ekf_list[-2]):
                        self.yaw_evaluated = np.arctan2(self.y_ekf_list[-1] - self.y_ekf_list[-2],
                                                        self.x_ekf_list[-1] - self.x_ekf_list[-2])

                        if self.yaw_evaluated_prev > 0 and self.yaw_evaluated <= -np.pi / 2:
                            self.yaw_evaluated += (2 * np.pi)
                            self.yaw_evaluated_list.append(self.yaw_evaluated)
                        else:
                            self.yaw_evaluated_list.append(self.yaw_evaluated)
                    else:
                        self.yaw_evaluated_list.append(self.yaw_evaluated_prev)

                    self.yaw_evaluated_prev = self.yaw_evaluated

                    # Gyros information
                    self.gyr_x = (float(self.line_list[67]) - self.init_gyr_x) / self.GYROSCOPE_SENSITIVITY
                    self.gyr_x_list.append(self.gyr_x)

                    self.gyr_y = (float(self.line_list[68]) - self.init_gyr_y) / self.GYROSCOPE_SENSITIVITY
                    self.gyr_y_list.append(self.gyr_y)

                    self.gyr_z = (float(self.line_list[69]) - self.init_gyr_z) / self.GYROSCOPE_SENSITIVITY
                    self.gyr_z_list.append(self.gyr_z)

                    self.acc_x = -float(self.line_list[64]) / self.ACCELEROMETER_SENSITIVITY
                    self.acc_x_list.append(self.acc_x)

                    self.acc_y = -float(self.line_list[65]) / self.ACCELEROMETER_SENSITIVITY
                    self.acc_y_list.append(self.acc_y)

                    self.acc_z = -float(self.line_list[66]) / self.ACCELEROMETER_SENSITIVITY
                    self.acc_z_list.append(self.acc_z)

                    at = int(self.line_list[0]) - self.prev_time

                    # Turning around the Y axis results in a vector on the X-axis
                    rollAcc = np.arctan2(self.acc_x, -self.acc_z) * 180.0 / np.pi;
                    self.roll += self.gyr_y * at / 1000000.0
                    self.roll = (self.roll * 0.98 + rollAcc * 0.02);
                    roll_corrected = -(self.roll + 4.15) * 6.5 / 15.0  # Roll compensated
                    self.roll_list.append(roll_corrected * np.pi / 180.0)  # Radians

                    pitchAcc = np.arctan2(self.acc_y, -self.acc_z) * 180.0 / np.pi;
                    self.pitch += self.gyr_x * at / 1000000.0
                    self.pitch = (self.pitch * 0.98 + pitchAcc * 0.02);
                    pitch_corrected = self.pitch + 0.04
                    self.pitch_list.append(pitch_corrected * np.pi / 180.0)  # Radians

                    data = {'timeStamp_ms': (self.time - self.init_time_stamp) / 1000.0,
                            'GPS_time_corrected_ms': self.gps_time_corrected_ms[-1], 'x': self.x_ekf_list[-1],
                            'y': self.y_ekf_list[-1], 'yaw': self.yaw_list_corrected[-1],
                            'pitch': (self.pitch_list[-1]) - 0.04,
                            'roll': self.roll_list[-1], 'calidad_gps': self.calidad_gps_list[-1]}

                    self.prev_time = int(self.line_list[0])

                    if self.verbose == True and self.init_time_stamp != 0:
                        # print data
                        # print at/1000000.0
                        print [self.time / 1000.0, int(self.line_list[1]), self.gps_time_corrected_ms[-1],
                               (time_aux % 50000) / 1000.0, int(self.line_list[1][6:8])]
                        print '----------------'
            elif self.gps_time_ms >= stop_limit_time:

                # self.roll_filtered = medfilt(self.roll_list, 125)
                self.time_diff = np.diff(self.gps_time_corrected_ms)
                self.time_diff = np.append(self.time_diff, self.time_diff[-1])
                print self.time_diff

                f_x = interpolate.interp1d(self.gps_time_corrected_ms, self.x_ekf_list)
                f_y = interpolate.interp1d(self.gps_time_corrected_ms, self.y_ekf_list)
                f_yaw = interpolate.interp1d(self.gps_time_corrected_ms, self.yaw_list_corrected)
                f_roll = interpolate.interp1d(self.gps_time_corrected_ms, self.roll_list)
                f_pitch = interpolate.interp1d(self.gps_time_corrected_ms, self.pitch_list)
                f_calidad = interpolate.interp1d(self.gps_time_corrected_ms, self.calidad_gps_list)

                self.time_inter = np.arange(self.gps_time_corrected_ms[0], self.gps_time_corrected_ms[-1], 0.2)

                self.x_inter = np.append([], f_x(self.time_inter))
                self.y_inter = np.append([], f_y(self.time_inter))
                self.yaw_inter = np.append([], f_yaw(self.time_inter))
                self.rol_inter = np.append([], f_roll(self.time_inter))
                self.pit_inter = np.append([], f_pitch(self.time_inter))
                self.cal_inter = np.append([], f_calidad(self.time_inter))

                for index in range(0, len(self.time_inter)):
                    self.data = {'GPS_time_corrected_ms': self.time_inter[index], 'x': self.x_inter[index],
                                 'y': self.y_inter[index],
                                 'yaw': self.yaw_inter[index], 'pitch': self.pit_inter[index],
                                 'roll': self.rol_inter[index],
                                 'quality': self.cal_inter[index]}

                    self.data_queue.append(self.data.copy())

                self.is_finished = True
                self.total_duration_gps = int(self.time_inter[-1]) - self.time_inter[0]
                print 'GPS Init time (us): ' + str(self.time_inter[0]) + ' GPS Final time (ms): ' + str(
                    self.time_inter[-1]) + ' Duration (minutes): ' + str(self.total_duration_gps / (60.0 * 1000.0))


                if self.show_data_flag == True:
                    self.show_data()

    def get_ekf_data(self):
        return self.data_queue


    def show_data(self):


        fig = plt.figure(1)

        ax = fig.gca(projection='3d')
        lines = ax.plot(self.x_list,self.y_list,self.gps_time_corrected_ms, picker=2.5)
        lines1 = ax.plot(self.x_ekf_list, self.y_list, self.gps_time_corrected_ms, color='red', picker=2.5)
        ax.set_xlabel("lon (m) x")
        ax.set_ylabel("lat (m) y")
        plt.title('Trajectory')
        plt.legend(['GPS', 'EKF'], loc='upper left')

        formatter = matplotlib.ticker.FuncFormatter(lambda ms, x: datetime.datetime.utcfromtimestamp(ms /1000.).strftime('%H:%M:%S'))
        ax.zaxis.set_major_formatter(formatter)

        fmt = matplotlib.ticker.FuncFormatter(lambda ms, x: datetime.datetime.utcfromtimestamp(ms /1000.).strftime('%H:%M:%S.%f')[:-3])
        datacursor(lines, formatter = lambda **kwargs: 'longitud(m): {x:.3f}\n'.format(**kwargs) +  'latitud(m): {y:.3f}\n'.format(**kwargs)+ '  hora: ' + fmt(kwargs.get('z')))
        datacursor(lines1, formatter = lambda **kwargs: 'longitud(m): {x:.3f}\n'.format(**kwargs) +  'latitud(m): {y:.3f}\n'.format(**kwargs)+ '  hora: ' + fmt(kwargs.get('z')))


        plt.figure(2)
        plt.subplot(221)
        plt.plot(self.gps_time_corrected_ms, self.roll_list)
        plt.plot(self.time_inter, self.rol_inter)
        plt.xlabel("time (s)")
        plt.ylabel("deg")
        plt.title('Roll')
        plt2 = plt.twinx()
        aux = np.array(self.roll_list, dtype=np.float32) * 180.0 / np.pi
        plt2.plot(self.gps_time_corrected_ms, aux, color='red')
        plt2.set_ylabel('deg', color='r')
        plt2.tick_params('y', colors='r')

        plt.subplot(222)
        plt.plot(self.time_inter, self.yaw_inter)
        plt.plot(self.gps_time_corrected_ms, self.yaw_list)
        plt.xlabel("time (s)")
        plt.ylabel("deg")
        plt.title('Yaw')

        plt.subplot(223)
        plt.plot(self.gps_time_corrected_ms, self.time_diff)
        plt.xlabel("time (s)")
        plt.ylabel("deg")
        plt.title('Gyro Y')

        plt.subplot(224)
        plt.plot(self.x_list, self.y_list)
        plt.plot(self.x_ekf_list, self.y_ekf_list, color='red')
        plt.plot(self.x_list[0], self.y_list[0], 'x', markersize=12)
        plt.plot(self.x_ekf_list[0], self.y_ekf_list[0], 'x', markersize=12)
        plt.plot(self.x_calidad_list, self.y_calidad_list, color='blue', linewidth=2)
        # plt.ylim([-250, 50])
        plt.xlabel("lon (m)")
        plt.ylabel("lat (m)")
        plt.title('Trajectory')
        plt.legend(['GPS', 'EKF'], loc='upper left')

        plt.tight_layout()
        plt.show()


import matplotlib.pyplot as plt
from scipy import interpolate



if __name__ == "__main__":

    # x = np.array([57967241.642], dtype=np.float32)
    # x = np.append(x, 57967251)
    # x = np.append(x, 57967261)
    # y = np.array([1.642, 1.8, 2], dtype=np.float32)
    # f = interpolate.interp1d(x, y)
    # xnew = np.arange(x[0], x[1], 1)
    # ynew = f(xnew)   # use interpolation function returned by `interp1d`
    # aux = xnew

    # xnew = np.arange(x[1], x[2], 1)
    # ynew = np.append(ynew, f(xnew))   # use interpolation function returned by `interp1d`
    # aux = np.append(aux, xnew)

    # print ynew

    # plt.plot(x, y, 'o', aux, ynew, 'x')
    # plt.show()
    # sys.exit()"""

    ekf = ekf_reader(True)
    ekf.read()



    if len(sys.argv) == 3:
        ekf = ekf_reader(sys.argv[1], sys.argv[2], True, True)
        ekf.read()
    elif len(sys.argv) == 2:
        ekf = ekf_reader(sys.argv[1], None, True, True)
        ekf.read()
    else:
        print 'Arguments are required'