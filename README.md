### EE152 Real-time Embedded Systems Final Project: Pan-Tompkins Algorithm for Real-Time ECG Analysis

#### Project Overview
This project implements and debugs a real-time signal-processing algorithm to detect QRS complexes in ECG waveforms using the Pan-Tompkins algorithm (PTC). Designed for an STM32 microcontroller with FreeRTOS, the system processes ECG signals, identifies QRS complexes, and calculates instantaneous heart rate. It includes debugging techniques to handle challenges in real-time biomedical signal processing.

#### Features
- **Hardware Integration**: Utilizes the STM32 Nucleo board, DACs, ADCs, and a seven-segment display for real-time ECG analysis.
- **Real-Time QRS Detection**: Identifies QRS complexes and signals them with audible beeps.
- **Heart Rate Monitoring**: Calculates and displays instantaneous heart rate on an LCD.
- **Debugging Tools**: Innovative use of DACs and GPIO pins for visualizing internal variables on an oscilloscope.
- **Live and Canned ECG Analysis**: Processes both prerecorded ECG data and live ECG signals.

#### Implementation Highlights
- **Debugging with Visualization**:
  - Outputs signals to DAC pins for real-time debugging with oscilloscopes.
  - Implements host-based debugging using a custom `analogRead` function and Python plotting tools.
- **Modular Codebase**: Includes modular FreeRTOS tasks and reusable libraries for ADC, DAC, GPIO, and UART operations.
- **Step-by-Step Progression**:
  - Initial setup with a canned ECG waveform.
  - Debugging PTC for various ECG datasets to handle algorithmic challenges.
  - Transition to live ECG data using an AD8232 preamplifier.

#### Future Directions
- Improve algorithm robustness to handle a wider variety of ECG patterns.
- Explore optimization techniques for real-time performance on embedded platforms.
- Extend visualization capabilities to include additional signals and metrics.

#### Repository Content
- **Source Code**: Modular C files for ECG signal processing and debugging. Includes a platformIO setup for use with a STM32L432KC nucleo board.
- **Data**: Various ECG data collected in 10 second intervals, trimmed down for useage with algorithm. 
- **Utilities**: Python scripts for signal plotting during host-based debugging under `Lab7_PTC.`

#### Authors
Gabby Beddard, Matt Bass (@mattjax16), Phaidra Martin (@pnmartin), Victoria Pontikes (@vpontikes)
