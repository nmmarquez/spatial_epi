# Disease Mapping  

### Why Disease mapping?

- Can be simple descriptive stats  
- Can have hypothesis  
- Background variability \(RR\)  
- Data comparability is difficult  
- Mortality is more sure  
- Morbidity is more interesting  

### Whats the scale 
 
- small scale means you risk bias and migration problems  
- large scale means you are smoothing over small scale variability  
- whatever the scale consider presenation 
    - cloropleth or isopleth  
    - color  
    - what are the cut points?  
    
### Smoothing models  

- estimates of SMR are often highly variable  
- how can we use all our geographical data to smooth out these estimates  
- the basic poisson model has no smoothing  
- RE effects model  
    - Poisson-Gamma  
    - Poisson-lognormal  
    - Poisson-lognormal-spatial
- all of these models take the form of  
$Y_{i} \theta_{i} \approx Poisson(E_{i} \theta_{i} )$  
$\theta_{i} = exp(\beta \times X) \delta_{i} \eta_{i}$
- sprinkle in covariates  
- what are the estimation techniques  
- remember that smoothing introduces bias  