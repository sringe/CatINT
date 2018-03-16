# `COMSOL`

There is a range of predefined COMSOL variables which will be written by default to the COMSOL input file. These can be used when defining new COMSOL variables or requesting specific output.

## General `COMSOL` Functions

There are some generalized `COMSOL` functions available for use, some of them listed here (`[var]` stands for the variable name):

|Functional Call | Equation |Meaning|
--- | --- | --- |
|`[var]x`, `d([var],x])` | $\frac{\mathrm{d}[\mathrm{var}]}{\mathrm{d}x}$|derivative with respect to x-direction
|`[var].at0(0,[var])`|$[\mathrm{var}](x=0)$|value of the variable at $x=0$
|`intop2([var]*(x<=dest(x)))`| $\displaystyle\int\limits_0^x [\mathrm{var}] \mathrm{d}x$ | integral of variable from 0 to x

`intop2` is defined as a domain integral as default integral.

## Species-Dependent Variables:

Species dependent variables can be requested by putting the species name `sp` into double brackets behind the physical property symbol:

|Variable Name `CatINT` | Variable Name `COMSOL` |Equation | Description |
 --- |--- |--- |--- |--- |--- |--- |--- |---
|`j[[sp]]`|`j1,j2,...` |$j_i(x)$|Flux of species `sp` |
|`c[[sp]]`|`cp1,cp2,...`|$c_i(x)$|Concentration of species `sp` |
|`ci[[sp]]`|`ci1,ci2,...`|$c_i(t=0)$|Initial Concentration


## Global Variables

|Variable Name `CatINT` | Variable Name `COMSOL` |Equation | Description |
 --- |--- |--- |--- |--- |--- |--- |--- |---
|`phi`|-- | Electrostatic Potential in Solution |
|`rho_s`|$\left(\Phi^\mathrm{M}-\Phi^\mathrm{PZC}\right)\cdot C_\mathrm{S}$  | Surface Charge Density |
|`rho_c`|$\rho_\mathrm{c}=F^2 \sum_i z_i^2 u_i c_i(x)$ | Electrolyte Conductivity|
|`i_el`|$i_\mathrm{el}=F \sum\limits_i z_i j_i$|Electrolyte Current Density

|`delta_phi_iR`|$\Delta \Phi_\mathrm{iR}(x)=\displaystyle\int\limits_0^x \frac{i(x)}{\rho_\mathrm{c}} \mathrm{d}x$|iR Drop in the Electrolyte as a function of distance to the electrode|

