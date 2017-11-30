import com.comsol.model.*;
import com.comsol.model.util.*;
/*
 *Transport Model as created by Transport Framework
 */
public class pnp_transport {
  public static Model run() {
    Model model = ModelUtil.create("Model");
    model.modelPath("/Users/sringe/software/catint/examples/CO_reduction");
    model.label("transport");
    model.comments("");
    model.param().set("D1", "5.273e-09[m^2/s]", "hydroxyl Diffusion coefficient");
    model.param().set("D2", "1.957e-09[m^2/s]", "potassium Diffusion coefficient");
    model.param().set("D3", "2.03e-09[m^2/s]", "carbon monoxide Diffusion coefficient");
    model.param().set("D4", "1.49e-09[m^2/s]", "C1-methane Diffusion coefficient");
    model.param().set("D5", "8.4e-10[m^2/s]", "C2-ethanol Diffusion coefficient");
    model.param().set("Z1", "-1", "hydroxyl ionic charge");
    model.param().set("Z2", "1", "potassium ionic charge");
    model.param().set("Z3", "0", "carbon monoxide ionic charge");
    model.param().set("Z4", "0", "C1-methane ionic charge");
    model.param().set("Z5", "0", "C2-ethanol ionic charge");
    model.param().set("k1", "0.95*A*exp(-max([Ga_CHO+0.5*V, Ga_CHOH+2*V])/RT)0.95**2*A*exp(-max([Ga_OCCOH+0.5*V, Ga_OCCO])/RT)[mol/m^2/s]", "hydroxyl flux");
    model.param().set("k2", "0.0[mol/m^2/s]", "potassium flux");
    model.param().set("k3", "0.95*A*exp(-max([Ga_CHO+0.5*V, Ga_CHOH+2*V])/RT)-2/12.0*0.95**2*A*exp(-max([Ga_OCCOH+0.5*V, Ga_OCCO])/RT)[mol/m^2/s]", "carbon monoxide flux");
    model.param().set("k4", "0.95*A*exp(-max([Ga_CHO+0.5*V, Ga_CHOH+2*V])/RT)[mol/m^2/s]", "C1-methane flux");
    model.param().set("k5", "0.95**2*A*exp(-max([Ga_OCCOH+0.5*V, Ga_OCCO])/RT)[mol/m^2/s]", "C2-ethanol flux");
    model.param().set("ci1", "1e-10[mol/m^3]", "hydroxyl bulk concentrations");
    model.param().set("ci2", "1e-10[mol/m^3]", "potassium bulk concentrations");
    model.param().set("ci3", "0.95[mol/m^3]", "carbon monoxide bulk concentrations");
    model.param().set("ci4", "0.0[mol/m^3]", "C1-methane bulk concentrations");
    model.param().set("ci5", "0.0[mol/m^3]", "C2-ethanol bulk concentrations");
    model.param().set("L_cell", "7.93e-05", "Cell length");
    model.param().set("epsilon", "lambdaD/L_cell", "Dimensionless Debye length scale");
    model.param().set("T", "298[K]", "Temperature");
    model.param().set("RT", "R_const*T", "Molar gas constant * Temperature");
    model.param().set("delta", "lambdaS/lambdaD", "Dimensionless Stern layer thickness");
    model.param().set("alphac", "0.5", "Cathodic charge transfer coefficient");
    model.param().set("alphaa", "1-alphac", "Anodic charge transfer coefficient");
    model.param().set("V", "0.0 [V]", "Potential");
    model.param().set("lambdaD", "0.000960882828277[m]", "Debye length");
    model.param().set("CS", "1.8e-05[F/cm^2]", "Stern layer capacitance");
    model.param().set("eps_r", "78.36", "relative permittivity");
    model.param().set("epsS", "epsilon0_const*2.0", "Stern layer effective permittivity");
    model.param().set("lambdaS", "epsS/CS", "Stern layer thickness");
    model.param().set("V", "0.2[V]", "Electrode potential");
    model.param().set("A", "1e+13",   "Exponential prefactor");
    model.param().set("Ga_CHOH", "2.37467774[eV]",   "CHOH Activation Energy");
    model.param().set("Ga_CHO", "1.11746219[eV]",   "CHO Activation Energy");
    model.param().set("Ga_OCCOH", "1.10495851[eV]",   "OCCOH Activation Energy");
    model.param().set("Ga_OCCO", "0.578959276[eV]",   "OCCO Activation Energy");
    model.param().set("flux_factor", "1", "factor scaling the flux");
    model.param().set("grid_factor", "40", "minimal number of x mesh points (boundary and domain)");
/*
 *MODEL DEFINITION
 */
    model.component().create("comp1", true);
    model.component("comp1").geom().create("geom1", 1);
    model.result().table().create("tbl1", "Table");
    model.component("comp1").mesh().create("mesh1");
    model.component("comp1").geom("geom1").create("i1", "Interval");
    model.component("comp1").geom("geom1").feature("i1").set("p2", "L_cell");
    model.component("comp1").geom("geom1").run();
/*
 *VARIABLES
 */
    model.component("comp1").variable().create("var1");
    model.component("comp1").variable("var1").set("deltaphi", "phiM-phi", "Metal - reaction plane potential difference");
    model.component("comp1").variable("var1").set("rho_s", "epsS*deltaphi/lambdaS", "Surface charge density");
    model.component("comp1").variable("var1").selection().geom("geom1", 0);
    model.component("comp1").variable("var1").selection().set(new int[]{1});
    model.component("comp1").variable().create("var2");
    model.component("comp1").variable("var2").set("phiM", "V", "Metal phase potential (ground)");
    model.component("comp1").variable("var2").selection().geom("geom1", 0);
    model.component("comp1").variable("var2").selection().set(new int[]{1});
    model.component("comp1").variable().create("var3");
    model.component("comp1").variable("var3").set("phiM", "0.0 [V]", "Metal phase potential (cell voltage)");
    model.component("comp1").variable("var3").selection().geom("geom1", 0);
    model.component("comp1").variable("var3").selection().set(new int[]{2});
/*
 *INTEGRATION
 */
    model.component("comp1").cpl().create("intop1", "Integration");
    model.component("comp1").cpl().create("intop2", "Integration");
    model.component("comp1").cpl("intop1").selection().geom("geom1", 0);
    model.component("comp1").cpl("intop1").selection().set(new int[]{2});
    model.component("comp1").cpl("intop2").selection().set(new int[]{1});
/*
 *ELECTROSTATICS
 */
    model.component("comp1").physics().create("es", "Electrostatics", "geom1");
    model.component("comp1").physics("es").field("electricpotential").field("phi");
    model.component("comp1").physics("es").create("sfcd1", "SurfaceChargeDensity", 0);
    model.component("comp1").physics("es").feature("sfcd1").selection().set(new int[]{1});
    model.component("comp1").physics("es").create("pot1", "ElectricPotential", 0);
    model.component("comp1").physics("es").feature("pot1").selection().set(new int[]{2});
    model.component("comp1").physics("es").create("pot2", "ElectricPotential", 0);
    model.component("comp1").physics("es").feature("pot2").selection().set(new int[]{1});
    model.component("comp1").physics("es").create("df1", "DisplacementField", 0);
    model.component("comp1").physics("es").feature("df1").selection().set(new int[]{1});
    model.component("comp1").physics().create("tds", "DilutedSpecies", "geom1");
    model.component("comp1").physics("tds").field("concentration").field("cp1");
    model.component("comp1").physics("tds").field("concentration").component(new String[]{"cp1", "cp2", "cp3", "cp4", "cp5"});
    model.component("comp1").physics("tds").create("fl1", "Fluxes", 0);
    model.component("comp1").physics("tds").feature("fl1").selection().set(new int[]{1});
    model.component("comp1").physics("tds").create("gconstr1", "GlobalConstraint", -1);
    model.component("comp1").physics("tds").create("conc1", "Concentration", 0);
    model.component("comp1").physics("tds").feature("conc1").selection().set(new int[]{2});
    model.component("comp1").physics("tds").create("reac1", "Reactions", 1);
    model.component("comp1").physics("tds").feature("reac1").selection().all();
    model.component("comp1").physics().create("ge", "GlobalEquations", "geom1");
    model.component("comp1").multiphysics().create("pc1", "PotentialCoupling", 1);
    model.component("comp1").multiphysics("pc1").selection().all();
    model.component("comp1").multiphysics().create("scdc1", "SpaceChargeDensityCoupling", 1);
    model.component("comp1").multiphysics("scdc1").selection().all();
    model.component("comp1").mesh("mesh1").create("edg1", "Edge");
    model.component("comp1").mesh("mesh1").feature("edg1").create("size1", "Size");
    model.component("comp1").mesh("mesh1").feature("edg1").create("size2", "Size");
    model.component("comp1").mesh("mesh1").feature("edg1").feature("size2").selection().geom("geom1", 0);
    model.component("comp1").mesh("mesh1").feature("edg1").feature("size2").selection().set(new int[]{1, 2});
    model.component("comp1").probe().create("pdom1", "DomainPoint");
    model.component("comp1").probe("pdom1").create("ppb2", "PointExpr");
    model.component("comp1").probe("pdom1").create("ppb3", "PointExpr");
    model.result().table("tbl1").label("Probe Table 1");
    model.component("comp1").variable("var3").active(false);
    model.component("comp1").view("view1").axis().set("xmin", 0);
    model.component("comp1").physics("es").feature("ccn1")
         .set("epsilonr", new String[][]{{"eps_r"}, {"0"}, {"0"}, {"0"}, {"eps_r"}, {"0"}, {"0"}, {"0"}, {"eps_r"}});
    model.component("comp1").physics("es").feature("sfcd1").set("rhoqs", "rho_s");
    model.component("comp1").physics("es").feature("pot2").set("V0", "V");
    model.component("comp1").physics("es").feature("pot2").active(false);
    model.component("comp1").physics("es").feature("df1").active(false);
    model.component("comp1").physics("es").feature("fp1").active(false);
    model.component("comp1").physics("tds").prop("ShapeProperty").set("order_concentration", 2);
    model.component("comp1").physics("tds").prop("TransportMechanism").set("Convection", false);
    model.component("comp1").physics("tds").prop("TransportMechanism").set("Migration", true);
/*
 *-> ADD diffusion coefficients
 */
    model.component("comp1").physics("tds").feature("cdm1").set("z", new String[][]{{"Z1"}, {"Z2"}, {"Z3"}, {"Z4"}, {"Z5"}});
    model.component("comp1").physics("tds").feature("cdm1").set("minput_temperature", "T");
    model.component("comp1").physics("tds").feature("cdm1").set("D_cp1", new String[][]{{"D1"}, {"0"}, {"0"}, {"0"}, {"D1"}, {"0"}, {"0"}, {"0"}, {"D1"}});
    model.component("comp1").physics("tds").feature("cdm1").set("D_cp2", new String[][]{{"D2"}, {"0"}, {"0"}, {"0"}, {"D2"}, {"0"}, {"0"}, {"0"}, {"D2"}});
    model.component("comp1").physics("tds").feature("cdm1").set("D_cp3", new String[][]{{"D3"}, {"0"}, {"0"}, {"0"}, {"D3"}, {"0"}, {"0"}, {"0"}, {"D3"}});
    model.component("comp1").physics("tds").feature("cdm1").set("D_cp4", new String[][]{{"D4"}, {"0"}, {"0"}, {"0"}, {"D4"}, {"0"}, {"0"}, {"0"}, {"D4"}});
    model.component("comp1").physics("tds").feature("cdm1").set("D_cp5", new String[][]{{"D5"}, {"0"}, {"0"}, {"0"}, {"D5"}, {"0"}, {"0"}, {"0"}, {"D5"}});
/*
 *-> ADD initial concentrations
 */
    model.component("comp1").physics("tds").feature("init1").set("initc", new String[][]{{"ci1"}, {"ci2"}, {"ci3"}, {"ci4"}, {"ci5"}});
/*
 *-> ADD BOUNDARIES
 */
    model.component("comp1").physics("tds").feature("fl1").set("species", new int[][]{{1}, {1}, {1}, {1}, {1}});
/*
 *-> ADD FLUXES
 */
    model.component("comp1").physics("tds").feature("fl1").set("N0", new String[][]{{"k0*flux_factor"}, {"k1*flux_factor"}, {"k2*flux_factor"}, {"k3*flux_factor"}, {"k4*flux_factor"}});
    model.component("comp1").physics("tds").feature("gconstr1").set("constraintExpression", "intop2(cm)-(cref*L)");
    model.component("comp1").physics("tds").feature("gconstr1").active(false);
    model.component("comp1").physics("tds").feature("conc1").set("species", new int[][]{{1}, {1}, {1}, {1}, {1}});
    model.component("comp1").physics("tds").feature("conc1").set("c0", new String[][]{{"ci1"}, {"ci2"}, {"ci3"}, {"ci4"}, {"ci5"}});
/*
 *REACTIONS
 */
    model.component("comp1").physics("ge").active(false);
    model.component("comp1").physics("ge").feature("ge1").set("name", "icell");
/*
 *GLOBAL EQUATIONS
 */
    model.component("comp1").physics("ge").feature("ge1").set("equation", "intop1(iloc)-icell");
    model.component("comp1").physics("ge").feature("ge1").set("description", "Cell Voltage");
    model.component("comp1").physics("ge").feature("ge1").set("DependentVariableQuantity", "electricpotential");
    model.component("comp1").physics("ge").feature("ge1").set("SourceTermQuantity", "currentdensity");
    model.component("comp1").physics("ge").feature("ge1")
         .label("Solve for potential that results in given current density icell");
    model.component("comp1").mesh("mesh1").feature("edg1").feature("size1").set("hauto", 1);
    model.component("comp1").mesh("mesh1").feature("edg1").feature("size1").set("custom", "on");
    model.component("comp1").mesh("mesh1").feature("edg1").feature("size1").set("hmax", "L_cell/grid_factor");
    model.component("comp1").mesh("mesh1").feature("edg1").feature("size1").set("hmaxactive", true);
    model.component("comp1").mesh("mesh1").feature("edg1").feature("size1").set("hnarrow", 5000);
    model.component("comp1").mesh("mesh1").feature("edg1").feature("size1").set("hnarrowactive", false);
    model.component("comp1").mesh("mesh1").feature("edg1").feature("size2").set("custom", "on");
    model.component("comp1").mesh("mesh1").feature("edg1").feature("size2").set("hmax", "lambdaD/grid_factor");
    model.component("comp1").mesh("mesh1").feature("edg1").feature("size2").set("hmaxactive", true);
    model.component("comp1").mesh("mesh1").feature("edg1").feature("size2").set("hmin", 3.26E-12);
    model.component("comp1").mesh("mesh1").feature("edg1").feature("size2").set("hminactive", false);
    model.component("comp1").mesh("mesh1").run();
/*
 *PROBES
 */
    model.component("comp1").probe("pdom1").set("bndsnap1", true);
    model.component("comp1").probe("pdom1").feature("ppb1").active(false);
    model.component("comp1").probe("pdom1").feature("ppb1").set("table", "tbl1");
    model.component("comp1").probe("pdom1").feature("ppb1").set("window", "window1");
    model.component("comp1").probe("pdom1").feature("ppb2").set("expr", "es.Jdx");
    model.component("comp1").probe("pdom1").feature("ppb2").set("unit", "A/m^2");
    model.component("comp1").probe("pdom1").feature("ppb2").set("descr", "Displacement current density, x component");
    model.component("comp1").probe("pdom1").feature("ppb2").set("table", "tbl1");
    model.component("comp1").probe("pdom1").feature("ppb2").set("window", "window1");
    model.component("comp1").probe("pdom1").feature("ppb3").active(false);
    model.component("comp1").probe("pdom1").feature("ppb3").set("expr", "es.Jdx-r*F_const*charge");
    model.component("comp1").probe("pdom1").feature("ppb3").set("unit", "A/m^2");
    model.component("comp1").probe("pdom1").feature("ppb3").set("descr", "es.Jdx-r*F_const*charge");
    model.component("comp1").probe("pdom1").feature("ppb3").set("table", "tbl1");
    model.component("comp1").probe("pdom1").feature("ppb3").set("window", "window1");
    model.component("comp1").physics("es").feature("ccn1").set("epsilonr_mat", "userdef");
/*
 *STUDIES
 */
    model.study().create("std0");
    model.study("std0").create("time", "Transient");
    model.study("std0").create("param", "Parametric");
    model.sol().create("sol1");
    model.sol("sol1").study("std0");
    model.sol("sol1").attach("std0");
    model.sol("sol1").create("st1", "StudyStep");
    model.sol("sol1").create("v1", "Variables");
    model.sol("sol0").create("t1", "Time");
    model.sol("sol1").feature("s1").create("fc1", "FullyCoupled");
    model.sol("sol1").feature("s1").create("d1", "Direct");
    model.sol("sol1").feature("s1").feature().remove("fcDef");
    model.study("std0").feature("time").set("tlist", "range(0.0,20.0,20001)");
    model.sol("sol1").attach("std0");
    model.sol("sol1").feature("v1").set("clist", new String[]{"0.0,20.0,20001"});
    model.sol("sol1").feature("t1").set("tlist", "0.0,20.0,20001");
    model.sol("sol1").feature("t1").set("maxorder", 2);
    model.sol("sol1").feature("t1").feature("fc1").set("maxiter", 8);
    model.sol("sol1").feature("t1").feature("fc1").set("damp", 0.9);
    model.sol("sol1").feature("t1").feature("fc1").set("jtech", "once");
    moidel.sol("sol1").runAll();
/*
 *OUTPUTS
 */
    model.result().export().create("data1", "Data");
    model.result().export("data1").set("expr", new String[]{"cp1", "cp2", "cp3", "cp4", "cp5"});
    model.result().export("data1").set("unit", new String[]{"mol/m^3", "mol/m^3", "mol/m^3", "mol/m^3", "mol/m^3"});
    model.result().export("data1").set("descr", new String[]{"Concentrations: hydroxyl, potassium, carbon monoxide, C1-methane, C2-ethanol"});
    model.result().export("data1").set("filename", "/Users/sringe/software/catint/examples/CO_reduction/results/concentrations.txt");
    model.result().export("data1").run();
    model.result().export().create("data2", "Data");
    model.result().export("data2").set("expr", new String[]{"phi", "es.Ex"});
    model.result().export("data2").set("unit", new String[]{"V", "V/m"});
    model.result().export("data2").set("descr", new String[]{"Potential, Field: phi, es.Ex"});
    model.result().export("data2").set("filename", "/Users/sringe/software/catint/examples/CO_reduction/results/electrostatics.txt");
    model.result().export("data2").run();
    return model;
  }

  public static Model run2(Model model) {
  return model;
  }

  public static void main(String[] args) {
    Model model = run();
    run2(model);
  }

}
