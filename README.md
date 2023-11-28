# SI100B project: Direction-Of-Arrival(DOA)

This project is about locating the sound source by a microphone array.

## Module 1: Narrow-band DoA estimation

Check narrowband_Project.m and simply click the run button.

It will read Observation_nb.mat and process it.

## Module 2: Broadband DoA estimation

Check the .m ducument which named 'Module2' and simply click the run button.

It will read and process Observation_wb.mat first and then Array_output_01~04.wav next.

If the running time of the program is a bit long, please wait patiently. It is not an endless loop or a bug.

## Module 3: Microphone array experiment

Check the .m ducument which named 'Module3' and simply check the run button.

It will read 音轨-1~4 and process it, which is consisted with two people's voices recorded with Aduacity. Actually the right channel order is 2143, which will be noticed in the code.

If the running time of the program is a bit long, please wait patiently. It is not an endless loop or a bug.

## BONUS PART

*bonus_video.mp4 is the demo of the bonus part, you can also see it in bilibili (https://www.bilibili.com/video/BV1c44y1E7re/).

Check the .m ducument which named 'bonus' and make sure you check the run button after connecting the device successfully.

You can control the number of sound source by adjusting the spinner component.

If matlab shows error at this stage, please close matlab and connect the microphone before opening matlab.

The name of the device we have appointed as '麦克风 (USB YDB01 Audio Effect)'.

After run, a gui is created, and after Start button is pushed, it will record 0.5s data and pause 0.15s to plot it in each cycle.

Some small error of estimation is normal.

If you want to stop recording, push the 'Stop' button.


P.S. If there is any problems not mentioned above(error, wrong answer,etc.), please contact us!
If you deleted any file by mistake, please unzip again or redownload the files.
