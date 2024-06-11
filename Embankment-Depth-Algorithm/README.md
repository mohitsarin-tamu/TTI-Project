# Calculating Embankment Depth

This folder contains the code used to find embankment depth. We have implemented the Peak Detection Algorithm to identify the peaks of the rock and then find the corresponding valleys, which are used to determine the embankment depth. Additionally, if we have the ground truth of the peaks, we can calculate the error in predicting the points.

Steps to run the code ('completecode.m'):

- Initially, if you want to run the complete code in one go, use 'completecode.m' and provide the correct folder name from where the data is to be read: data_csv = csvread('/Users/mohitsarin/Desktop/TTI/a240402b.csv').
- If you want to calculate the error and the ground truth file is available, update the ground truth file location in the code: data_excel = readtable('/Users/mohitsarin/Desktop/TTI/ground_truth.xlsx').
- Sample data is also attached herewith.

Steps to run the code divided into functions for clarity (code folder):
- This folder consists of 3 files, namely: main.m, peaks2.m, find_valleys.m
- Run the main.m file after setting the required paths such as csv_file which is the path for the data collected using laser and ground_truth which is the path of groud truth file.
- If the groud truths are not available leave that as empty. 
- Define mat_thickness in embankment depth function (default 0.5).
- Define x_length, and y_length as width of x and y axis as per the LASER used and the data gathered (default 4,2 respectively).
- Define inputs to peaks2 function 'MinPeakHeight' as minpeakheight, 'Threshold' as peaks2threshold, 'MinPeakDistance'as minpeakdistance (default 0.1,0.001,0.2 respectively)
- Define threshold_radius as the radius around which valley point is searched for (default 1.5).
- Define threshold to find the nearby predicted points from the ground truth (default 3.0).


