# ESSBAS
Extremely Simple SBAS 

Simple, fast and not very sophisticated. (You might also say not very good, but I find it useful to compare its results to those from fancier codes.) Some fun facts:
<ul>
  <li>It doesn't have temporal smoothing, beyond a minimal SVD truncation, so the results are kinda noisy</li>
  <li>It doesn't have atmospheric corrections (yet!) either</li>
  <li>It <i>does</i> have a phase closure test, but throws out any pixel that doesn't close from the whole solution</li>
  <li>I wrote the basics of it in an afternoon</li>
  <li>It is written in MATLAB, and requires the Mapping Toolbox</li>
  <li>It can read in interferograms from LICS, ARIA and GMTSAR</li>
  <li><b>I haven't uploaded all of it yet</b>. Soz.</li>
</ul>

It uses the functions <a href="https://www.mathworks.com/matlabcentral/fileexchange/25683-grdread2">grdread2</a> and <a href="https://www.mathworks.com/matlabcentral/fileexchange/26290-grdwrite2">grdwrite2</a> which one day I might be motivated to reverse engineer, but for now can be found on the <a href="https://www.mathworks.com/matlabcentral/fileexchange/">MATLAB Central File Exchange</a>. 
