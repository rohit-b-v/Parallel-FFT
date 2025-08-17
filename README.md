# 📌 Parallel-FFT

FFT is an improved algorithm for evalulating Fourier transform of signals, which otherwise had higher time complexity due to the number of heavy computations involved in it. FFT algorithm can be further sped up by integrating Parallel Computing tools into it. For the demonstration of the Parallel FFT, the well known Cooley and Tukey algorithm is used on the input data array and implemented as a parallel non-recursive version of FFT using OPENMP and MPI.

## ✨ Why This Project?
This project was built to **demonstrate practical skills** in:
- 🖥️ Programming: C++  
- 🛠️ Tools: OpenMP, MPI

## Screen shots
<img width="519" height="323" alt="image" src="https://github.com/user-attachments/assets/e86c61fb-2fc4-4b08-ae68-e4d7d928e05a" />
<img width="512" height="324" alt="image" src="https://github.com/user-attachments/assets/ad11b4e2-b913-4b70-8f2e-218c9e274748" />



## 🛠️ Tech Stack
- **Languages:** C++  
- **Other Tools:** MPI, OpenMP  

---

## 🚀 Features
- On running Cooley-Tukey algorithm on OpenMP and MPI there is an improvement in the execution time by 4 times when compared to serial execution.
- OpenMP was more efficient in parallelising the algorithm than MPI.
- There was approximately 2 times decrease in execution time when using OpenMP in comparison with MPI. 
---
