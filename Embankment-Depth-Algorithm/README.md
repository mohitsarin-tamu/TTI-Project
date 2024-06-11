# Calculating Embankment Depth

This folder contains the code used to find embankment depth. We have implemented the Peak Detection Algorithm to identify the peaks of the rock and then find the corresponding valleys, which are used to determine the embankment depth. Additionally, if we have the ground truth of the peaks, we can calculate the error in predicting the points.

Steps to run the code ('completecode.m'):

- Initially, if you want to run the complete code in one go, use 'completecode.m' and provide the correct folder name from where the data is to be read: data_csv = csvread('/Users/mohitsarin/Desktop/TTI/a240402b.csv').
- If you want to calculate the error and the ground truth file is available, update the ground truth file location in the code: data_excel = readtable('/Users/mohitsarin/Desktop/TTI/ground_truth.xlsx').
- Sample data is also attached herewith.

Steps to run the code divided into functions for clarity (code folder):

-
