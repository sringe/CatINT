/*
 * diffuse_double_layer_with_charge_transfer_nonstatic_2.java
 */

import com.comsol.model.*;
import com.comsol.model.util.*;

/** Model exported on Oct 31 2017, 11:22 by COMSOL 5.3.0.316. */
public class diffuse_double_layer_with_charge_transfer_nonstatic_2 {

  public static Model run() {
    Model model = ModelUtil.create("Model");

    model.modelPath("/Users/sringe/software/comsol/comsol_files/CarlosModel");

    model.label("diffuse_double_layer_with_charge_transfer_nonstatic_2.mph");

    model.comments("Diffuse Double Layer With Charge Transfer\n\nThis example is an extension of the Diffuse Double Layer model and demonstrates how incorporate charge transfer in a diffuse double layer described by the coupled Nernst-Planck-Poisson set of equations.");

    model.param().set("D1", "1.13990561575e-09[m^2/s]", "propol Diffusion coefficient");
    model.param().set("D2", "4.62363254756e-09[m^2/s]", "OH- Diffusion coefficient");
    model.param().set("D3", "7.36554397866e-10[m^2/s]", "metol Diffusion coefficient");
    model.param().set("D4", "1.67478440467e-09[m^2/s]", "CO2 Diffusion coefficient");
    model.param().set("D5", "1.78000646151e-09[m^2/s]", "CO Diffusion coefficient");
    model.param().set("D6", "9.66289221962e-10[m^2/s]", "etgly Diffusion coefficient");
    model.param().set("D7", "9.64535521015e-10[m^2/s]", "allyl Diffusion coefficient");
    model.param().set("D8", "3.94582713143e-09[m^2/s]", "H2 Diffusion coefficient");
    model.param().set("D9", "0.0[m^2/s]", "unknown Diffusion coefficient");
    model.param().set("D10", "1.71599637693e-09[m^2/s]", "K Diffusion coefficient");
    model.param().set("D11", "7.36554397866e-10[m^2/s]", "etol Diffusion coefficient");
    model.param().set("D12", "1.03906781128e-09[m^2/s]", "HCO3- Diffusion coefficient");
    model.param().set("D13", "1.30650720574e-09[m^2/s]", "CH4 Diffusion coefficient");
    model.param().set("D14", "1.63971038573e-09[m^2/s]", "C2H4 Diffusion coefficient");
    model.param().set("D15", "1.27494058869e-09[m^2/s]", "HCOO- Diffusion coefficient");
    model.param().set("D16", "8.09332987179e-10[m^2/s]", "CO32- Diffusion coefficient");
    model.param().set("D17", "9.54890165805e-10[m^2/s]", "acet Diffusion coefficient");
    model.param().set("Z1", "0", "propol charge");
    model.param().set("Z2", "-1", "OH- charge");
    model.param().set("Z3", "0", "metol charge");
    model.param().set("Z4", "0", "CO2 charge");
    model.param().set("Z5", "0", "CO charge");
    model.param().set("Z6", "0", "etgly charge");
    model.param().set("Z7", "0", "allyl charge");
    model.param().set("Z8", "0", "H2 charge");
    model.param().set("Z9", "0", "unknown charge");
    model.param().set("Z10", "1", "K charge");
    model.param().set("Z11", "0", "etol charge");
    model.param().set("Z12", "-1", "HCO3- charge");
    model.param().set("Z13", "0", "CH4 charge");
    model.param().set("Z14", "0", "C2H4 charge");
    model.param().set("Z15", "-1", "HCOO- charge");
    model.param().set("Z16", "-2", "CO32- charge");
    model.param().set("Z17", "-1", "acet charge");
    model.param().set("k1", "1.32790476537e-07[mol/m^2/s]", "propol flux");
    model.param().set("k2", "0.00187962930218[mol/m^2/s]", "OH- flux");
    model.param().set("k3", "3.38241305932e-08[mol/m^2/s]", "metol flux");
    model.param().set("k4", "-0.000132386711388[mol/m^2/s]", "CO2 flux");
    model.param().set("k5", "2.60588829897e-06[mol/m^2/s]", "CO flux");
    model.param().set("k6", "2.39269069282e-08[mol/m^2/s]", "etgly flux");
    model.param().set("k7", "0.0[mol/m^2/s]", "allyl flux");
    model.param().set("k8", "0.000464243493372[mol/m^2/s]", "H2 flux");
    model.param().set("k9", "0.0[mol/m^2/s]", "unknown flux");
    model.param().set("k10", "0.0[mol/m^2/s]", "K flux");
    model.param().set("k11", "2.50704426004e-06[mol/m^2/s]", "etol flux");
    model.param().set("k12", "0.0[mol/m^2/s]", "HCO3- flux");
    model.param().set("k13", "9.21248907348e-05[mol/m^2/s]", "CH4 flux");
    model.param().set("k14", "1.39220331191e-05[mol/m^2/s]", "C2H4 flux");
    model.param().set("k15", "4.16043545663e-06[mol/m^2/s]", "HCOO- flux");
    model.param().set("k16", "0.0[mol/m^2/s]", "CO32- flux");
    model.param().set("k17", "7.86463827995e-08[mol/m^2/s]", "acet flux");
    model.param().set("ci1", "0.0[mol/m^3]", "propol bulk concentrations");
    model.param().set("ci2", "6.58543643499e-05[mol/m^3]", "OH- bulk concentrations");
    model.param().set("ci3", "0.0[mol/m^3]", "metol bulk concentrations");
    model.param().set("ci4", "34.19[mol/m^3]", "CO2 bulk concentrations");
    model.param().set("ci5", "0.0[mol/m^3]", "CO bulk concentrations");
    model.param().set("ci6", "0.0[mol/m^3]", "etgly bulk concentrations");
    model.param().set("ci7", "0.0[mol/m^3]", "allyl bulk concentrations");
    model.param().set("ci8", "0.0[mol/m^3]", "H2 bulk concentrations");
    model.param().set("ci9", "0.0[mol/m^3]", "unknown bulk concentrations");
    model.param().set("ci10", "100.0[mol/m^3]", "K bulk concentrations");
    model.param().set("ci11", "0.0[mol/m^3]", "etol bulk concentrations");
    model.param().set("ci12", "99.9692958402[mol/m^3]", "HCO3- bulk concentrations");
    model.param().set("ci13", "0.0[mol/m^3]", "CH4 bulk concentrations");
    model.param().set("ci14", "0.0[mol/m^3]", "C2H4 bulk concentrations");
    model.param().set("ci15", "0.0[mol/m^3]", "HCOO- bulk concentrations");
    model.param().set("ci16", "0.0307041597553[mol/m^3]", "CO32- bulk concentrations");
    model.param().set("ci17", "0.0[mol/m^3]", "acet bulk concentrations");
    model.param().set("L_cell", "7.93e-05", "Cell length");
    model.param().set("epsilon", "lambdaD/L_cell", "Dimensionless Debye length scale");
    model.param().set("T", "298[K]", "Temperature");
    model.param().set("RT", "R_const*T", "Molar gas constant * Temperature");
    model.param()
         .set("delta", "lambdaS/lambdaD", "Dimensionless Stern layer thicknesscM 1[mol/m^3] Metal reference concentration");
    model.param().set("alphac", "0.5", "Cathodic charge transfer coefficient");
    model.param().set("alphaa", "1-alphac", "Anodic charge transfer coefficient");
    model.param().set("V", "-0.2 [V]", "Potential");
    model.param().set("lambdaD", "9.60661473344e-10[m]", "Debye length");
    model.param().set("CS", "18*1e-6[F/cm^2]", "Stern layer capacitance");
    model.param().set("eps_r", "78.36", "relative permittivity");
    model.param().set("epsS", "epsilon0_const*2", "Stern layer effective permittivity");
    model.param().set("lambdaS", "epsS/CS", "Stern layer thicknessV 0.2[V] Electrode potential");
    model.param().set("k1f", "1/1000*(5.93e3) [1/s*m^3/mol]", "rate constant: CO2 + OH- -> HCO3-");
    model.param().set("k1r", "1/(4.44e7)*(5.93e3) [1/s]", "rate constant: HCO3- -> CO2 + OH-");
    model.param().set("k2f", "1/1000.0*(1.0e8) [1/s*m^3/mol]", "rate constant: HCO3- + OH- -> CO32- + H2O");
    model.param().set("k2r", "1/(4.66e3)*(1.0e8) [1/s*m^3/mol]", "rate constant: CO32- + H2O -> HCO3- + OH-");
    model.param().set("flux_factor", "1", "factor scaling the flux");
    model.param().set("grid_factor", "40", "factor increasing the grid density from its default");

    model.component().create("comp1", true);

    model.component("comp1").geom().create("geom1", 1);

    model.result().table().create("tbl1", "Table");

    model.component("comp1").mesh().create("mesh1");

    model.component("comp1").geom("geom1").create("i1", "Interval");
    model.component("comp1").geom("geom1").feature("i1").set("p2", "L_cell");
    model.component("comp1").geom("geom1").run();

    model.component("comp1").variable().create("var1");
    model.component("comp1").variable("var1")
         .set("deltaphi", "phiM-phi", "Metal - reaction plane potential difference");
    model.component("comp1").variable("var1").set("rho_s", "epsS*deltaphi/lambdaS", "Surface charge density");
    model.component("comp1").variable("var1").selection().geom("geom1", 0);
    model.component("comp1").variable("var1").selection().set(new int[]{1});
    model.component("comp1").variable().create("var2");
    model.component("comp1").variable("var2").set("phiM", "V", "Metal phase potential (ground)");
    model.component("comp1").variable("var2").selection().geom("geom1", 0);
    model.component("comp1").variable("var2").selection().set(new int[]{1});
    model.component("comp1").variable().create("var3");
    model.component("comp1").variable("var3").set("phiM", "1 [V]", "Metal phase potential (cell voltage)");
    model.component("comp1").variable("var3").selection().geom("geom1", 0);
    model.component("comp1").variable("var3").selection().set(new int[]{2});

    model.component("comp1").cpl().create("intop1", "Integration");
    model.component("comp1").cpl().create("intop2", "Integration");
    model.component("comp1").cpl("intop1").selection().geom("geom1", 0);
    model.component("comp1").cpl("intop1").selection().set(new int[]{2});
    model.component("comp1").cpl("intop2").selection().set(new int[]{1});

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
    model.component("comp1").physics("es").create("fp1", "FloatingPotential", 0);
    model.component("comp1").physics().create("tds", "DilutedSpecies", "geom1");
    model.component("comp1").physics("tds").field("concentration").field("cp1");
    model.component("comp1").physics("tds").field("concentration")
         .component(new String[]{"cp1", "cp2", "cp3", "cp4", "cp5", "cp6", "cp7", "cp8", "cp9", "cp10", 
         "cp11", "cp12", "cp13", "cp14", "cp15", "cp16", "cp17"});
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
    model.component("comp1").physics("tds").feature("cdm1")
         .set("z", new String[][]{{"Z1"}, 
         {"Z2"}, 
         {"Z3"}, 
         {"Z4"}, 
         {"Z5"}, 
         {"Z6"}, 
         {"Z7"}, 
         {"Z8"}, 
         {"Z9"}, 
         {"Z10"}, 
         {"Z11"}, 
         {"Z12"}, 
         {"Z13"}, 
         {"Z14"}, 
         {"Z15"}, 
         {"Z16"}, 
         {"Z17"}});
    model.component("comp1").physics("tds").feature("cdm1").set("minput_temperature", "T");
    model.component("comp1").physics("tds").feature("cdm1")
         .set("D_cp5", new String[][]{{"D5"}, {"0"}, {"0"}, {"0"}, {"D5"}, {"0"}, {"0"}, {"0"}, {"D5"}});
    model.component("comp1").physics("tds").feature("cdm1")
         .set("D_cp6", new String[][]{{"D6"}, {"0"}, {"0"}, {"0"}, {"D6"}, {"0"}, {"0"}, {"0"}, {"D6"}});
    model.component("comp1").physics("tds").feature("cdm1")
         .set("D_cp7", new String[][]{{"D7"}, {"0"}, {"0"}, {"0"}, {"D7"}, {"0"}, {"0"}, {"0"}, {"D7"}});
    model.component("comp1").physics("tds").feature("cdm1")
         .set("D_cp8", new String[][]{{"D8"}, {"0"}, {"0"}, {"0"}, {"D8"}, {"0"}, {"0"}, {"0"}, {"D8"}});
    model.component("comp1").physics("tds").feature("cdm1")
         .set("D_cp9", new String[][]{{"D9"}, {"0"}, {"0"}, {"0"}, {"D9"}, {"0"}, {"0"}, {"0"}, {"D9"}});
    model.component("comp1").physics("tds").feature("cdm1")
         .set("D_cp10", new String[][]{{"D10"}, {"0"}, {"0"}, {"0"}, {"D10"}, {"0"}, {"0"}, {"0"}, {"D10"}});
    model.component("comp1").physics("tds").feature("cdm1")
         .set("D_cp11", new String[][]{{"D11"}, {"0"}, {"0"}, {"0"}, {"D11"}, {"0"}, {"0"}, {"0"}, {"D11"}});
    model.component("comp1").physics("tds").feature("cdm1")
         .set("D_cp12", new String[][]{{"D12"}, {"0"}, {"0"}, {"0"}, {"D12"}, {"0"}, {"0"}, {"0"}, {"D12"}});
    model.component("comp1").physics("tds").feature("cdm1")
         .set("D_cp13", new String[][]{{"D13"}, {"0"}, {"0"}, {"0"}, {"D13"}, {"0"}, {"0"}, {"0"}, {"D13"}});
    model.component("comp1").physics("tds").feature("cdm1")
         .set("D_cp14", new String[][]{{"D14"}, {"0"}, {"0"}, {"0"}, {"D14"}, {"0"}, {"0"}, {"0"}, {"D14"}});
    model.component("comp1").physics("tds").feature("cdm1")
         .set("D_cp15", new String[][]{{"D15"}, {"0"}, {"0"}, {"0"}, {"D15"}, {"0"}, {"0"}, {"0"}, {"D15"}});
    model.component("comp1").physics("tds").feature("cdm1")
         .set("D_cp16", new String[][]{{"D16"}, {"0"}, {"0"}, {"0"}, {"D16"}, {"0"}, {"0"}, {"0"}, {"D16"}});
    model.component("comp1").physics("tds").feature("cdm1")
         .set("D_cp17", new String[][]{{"D17"}, {"0"}, {"0"}, {"0"}, {"D17"}, {"0"}, {"0"}, {"0"}, {"D17"}});
    model.component("comp1").physics("tds").feature("cdm1")
         .set("D_cp1", new String[][]{{"D1"}, {"0"}, {"0"}, {"0"}, {"D1"}, {"0"}, {"0"}, {"0"}, {"D1"}});
    model.component("comp1").physics("tds").feature("cdm1")
         .set("D_cp2", new String[][]{{"D2"}, {"0"}, {"0"}, {"0"}, {"D2"}, {"0"}, {"0"}, {"0"}, {"D2"}});
    model.component("comp1").physics("tds").feature("cdm1")
         .set("D_cp3", new String[][]{{"D3"}, {"0"}, {"0"}, {"0"}, {"D3"}, {"0"}, {"0"}, {"0"}, {"D3"}});
    model.component("comp1").physics("tds").feature("cdm1")
         .set("D_cp4", new String[][]{{"D4"}, {"0"}, {"0"}, {"0"}, {"D4"}, {"0"}, {"0"}, {"0"}, {"D4"}});
    model.component("comp1").physics("tds").feature("init1")
         .set("initc", new String[][]{{"ci1"}, 
         {"ci2"}, 
         {"ci3"}, 
         {"ci4"}, 
         {"ci5"}, 
         {"ci6"}, 
         {"ci7"}, 
         {"ci8"}, 
         {"ci9"}, 
         {"ci10"}, 
         {"ci11"}, 
         {"ci12"}, 
         {"ci13"}, 
         {"ci14"}, 
         {"ci15"}, 
         {"ci16"}, 
         {"ci17"}});
    model.component("comp1").physics("tds").feature("fl1")
         .set("species", new int[][]{{1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}});
    model.component("comp1").physics("tds").feature("fl1")
         .set("N0", new String[][]{{"k1*flux_factor"}, 
         {"k2*flux_factor"}, 
         {"k3*flux_factor"}, 
         {"k4*flux_factor"}, 
         {"k5*flux_factor"}, 
         {"k6*flux_factor"}, 
         {"k7*flux_factor"}, 
         {"k8*flux_factor"}, 
         {"k9*flux_factor"}, 
         {"k10*flux_factor"}, 
         {"k11*flux_factor"}, 
         {"k12*flux_factor"}, 
         {"k13*flux_factor"}, 
         {"k14*flux_factor"}, 
         {"k15*flux_factor"}, 
         {"k16*flux_factor"}, 
         {"k17*flux_factor"}});
    model.component("comp1").physics("tds").feature("gconstr1").set("constraintExpression", "intop2(cm)-(cref*L)");
    model.component("comp1").physics("tds").feature("gconstr1").active(false);
    model.component("comp1").physics("tds").feature("conc1")
         .set("species", new int[][]{{1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}});
    model.component("comp1").physics("tds").feature("conc1")
         .set("c0", new String[][]{{"ci1"}, 
         {"ci2"}, 
         {"ci3"}, 
         {"ci4"}, 
         {"ci5"}, 
         {"ci6"}, 
         {"ci7"}, 
         {"ci8"}, 
         {"ci9"}, 
         {"ci10"}, 
         {"ci11"}, 
         {"ci12"}, 
         {"ci13"}, 
         {"ci14"}, 
         {"ci15"}, 
         {"ci16"}, 
         {"ci17"}});
    model.component("comp1").physics("tds").feature("reac1")
         .set("R_cp2", "-k1f*cp4*cp2+k1r*cp12-k2f*cp12*cp2+k2r*cp16");
    model.component("comp1").physics("tds").feature("reac1").set("R_cp4", "-k1f*cp4*cp2+k1r*cp12");
    model.component("comp1").physics("tds").feature("reac1")
         .set("R_cp12", "k1f*cp4*cp2-k1r*cp12-k2f*cp12*cp2+k2r*cp16");
    model.component("comp1").physics("tds").feature("reac1").set("R_cp16", "k2f*cp12*cp2-k2r*cp16");
    model.component("comp1").physics("tds").feature("reac1").active(false);
    model.component("comp1").physics("ge").active(false);
    model.component("comp1").physics("ge").feature("ge1").set("name", "icell");
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

    model.study().create("std1");
    model.study("std1").create("stat", "Stationary");
    model.study().create("std2");
    model.study("std2").create("param", "Parametric");
    model.study("std2").create("time", "Transient");

    model.sol().create("sol1");
    model.sol("sol1").study("std1");
    model.sol("sol1").attach("std1");
    model.sol("sol1").create("st1", "StudyStep");
    model.sol("sol1").create("v1", "Variables");
    model.sol("sol1").create("s1", "Stationary");
    model.sol("sol1").feature("s1").create("fc1", "FullyCoupled");
    model.sol("sol1").feature("s1").create("d1", "Direct");
    model.sol("sol1").feature("s1").feature().remove("fcDef");
    model.sol().create("sol2");
    model.sol("sol2").study("std2");
    model.sol("sol2").attach("std2");
    model.sol("sol2").create("st1", "StudyStep");
    model.sol("sol2").create("v1", "Variables");
    model.sol("sol2").create("t1", "Time");
    model.sol("sol2").feature("t1").create("fc1", "FullyCoupled");
    model.sol("sol2").feature("t1").create("d1", "Direct");
    model.sol("sol2").feature("t1").feature().remove("fcDef");
    model.sol().create("sol5");
    model.sol("sol5").study("std2");
    model.sol("sol5").label("Parametric Solutions 3");

    model.batch().create("p1", "Parametric");
    model.batch("p1").create("so1", "Solutionseq");
    model.batch("p1").study("std2");

    model.result().dataset().create("cpt1", "CutPoint1D");
    model.result().dataset().create("dset6", "Solution");
    model.result().dataset("dset2").set("probetag", "pdom1");
    model.result().dataset("cpt1").set("probetag", "pdom1");
    model.result().dataset("cpt1").set("data", "dset2");
    model.result().dataset("dset3").set("solution", "sol2");
    model.result().dataset("dset6").set("solution", "sol5");
    model.result().numerical().create("pev1", "EvalPoint");
    model.result().numerical().create("pev2", "EvalPoint");
    model.result().numerical().create("pev3", "EvalPoint");
    model.result().numerical("pev1").set("probetag", "pdom1/ppb1");
    model.result().numerical("pev2").set("probetag", "pdom1/ppb2");
    model.result().numerical("pev3").set("probetag", "pdom1/ppb3");
    model.result().create("pg8", "PlotGroup1D");
    model.result().create("pg9", "PlotGroup1D");
    model.result().create("pg10", "PlotGroup1D");
    model.result().create("pg11", "PlotGroup1D");
    model.result().create("pg12", "PlotGroup1D");
    model.result().create("pg13", "PlotGroup1D");
    model.result().create("pg14", "PlotGroup1D");
    model.result().create("pg15", "PlotGroup1D");
    model.result().create("pg16", "PlotGroup1D");
    model.result("pg8").set("probetag", "window1_default");
    model.result("pg8").create("tblp1", "Table");
    model.result("pg8").feature("tblp1").set("probetag", "pdom1/ppb1,pdom1/ppb2,pdom1/ppb3");
    model.result("pg9").set("data", "dset6");
    model.result("pg9").create("lngr1", "LineGraph");
    model.result("pg9").create("lngr2", "LineGraph");
    model.result("pg9").create("lngr3", "LineGraph");
    model.result("pg9").feature("lngr1").set("xdata", "expr");
    model.result("pg9").feature("lngr1").selection().set(new int[]{1});
    model.result("pg9").feature("lngr2").set("xdata", "expr");
    model.result("pg9").feature("lngr2").selection().all();
    model.result("pg10").set("data", "dset3");
    model.result("pg10").create("lngr1", "LineGraph");
    model.result("pg10").feature("lngr1").selection().set(new int[]{1});
    model.result("pg11").set("data", "dset3");
    model.result("pg11").create("lngr1", "LineGraph");
    model.result("pg11").create("lngr2", "LineGraph");
    model.result("pg11").create("lngr3", "LineGraph");
    model.result("pg11").create("lngr4", "LineGraph");
    model.result("pg11").create("lngr5", "LineGraph");
    model.result("pg11").create("lngr6", "LineGraph");
    model.result("pg11").create("lngr7", "LineGraph");
    model.result("pg11").create("lngr8", "LineGraph");
    model.result("pg11").create("lngr9", "LineGraph");
    model.result("pg11").create("lngr10", "LineGraph");
    model.result("pg11").create("lngr11", "LineGraph");
    model.result("pg11").feature("lngr1").set("xdata", "expr");
    model.result("pg11").feature("lngr1").selection().set(new int[]{1});
    model.result("pg11").feature("lngr2").set("xdata", "expr");
    model.result("pg11").feature("lngr2").selection().set(new int[]{1});
    model.result("pg11").feature("lngr3").set("xdata", "expr");
    model.result("pg11").feature("lngr3").selection().set(new int[]{1});
    model.result("pg11").feature("lngr4").set("xdata", "expr");
    model.result("pg11").feature("lngr4").selection().set(new int[]{1});
    model.result("pg11").feature("lngr5").set("xdata", "expr");
    model.result("pg11").feature("lngr5").selection().set(new int[]{1});
    model.result("pg11").feature("lngr6").set("xdata", "expr");
    model.result("pg11").feature("lngr6").selection().set(new int[]{1});
    model.result("pg11").feature("lngr7").set("xdata", "expr");
    model.result("pg11").feature("lngr7").selection().set(new int[]{1});
    model.result("pg11").feature("lngr8").set("xdata", "expr");
    model.result("pg11").feature("lngr8").selection().set(new int[]{1});
    model.result("pg11").feature("lngr9").set("xdata", "expr");
    model.result("pg11").feature("lngr9").selection().set(new int[]{1});
    model.result("pg11").feature("lngr10").set("xdata", "expr");
    model.result("pg11").feature("lngr10").selection().set(new int[]{1});
    model.result("pg11").feature("lngr11").set("xdata", "expr");
    model.result("pg11").feature("lngr11").selection().set(new int[]{1});
    model.result("pg12").set("data", "dset6");
    model.result("pg12").create("lngr1", "LineGraph");
    model.result("pg12").create("lngr2", "LineGraph");
    model.result("pg12").create("lngr3", "LineGraph");
    model.result("pg12").create("lngr5", "LineGraph");
    model.result("pg12").create("lngr6", "LineGraph");
    model.result("pg12").feature("lngr1").set("xdata", "expr");
    model.result("pg12").feature("lngr1").selection().set(new int[]{1});
    model.result("pg12").feature("lngr2").set("xdata", "expr");
    model.result("pg12").feature("lngr2").selection().set(new int[]{1});
    model.result("pg12").feature("lngr3").set("xdata", "expr");
    model.result("pg12").feature("lngr3").selection().set(new int[]{1});
    model.result("pg12").feature("lngr5").set("xdata", "expr");
    model.result("pg12").feature("lngr5").selection().set(new int[]{1});
    model.result("pg12").feature("lngr6").set("xdata", "expr");
    model.result("pg12").feature("lngr6").selection().set(new int[]{1});
    model.result("pg13").set("data", "dset6");
    model.result("pg13").create("lngr1", "LineGraph");
    model.result("pg13").create("lngr2", "LineGraph");
    model.result("pg13").create("lngr3", "LineGraph");
    model.result("pg13").feature("lngr1").set("xdata", "expr");
    model.result("pg13").feature("lngr1").selection().set(new int[]{1});
    model.result("pg13").feature("lngr2").set("xdata", "expr");
    model.result("pg13").feature("lngr2").selection().set(new int[]{1});
    model.result("pg13").feature("lngr3").set("xdata", "expr");
    model.result("pg13").feature("lngr3").selection().set(new int[]{1});
    model.result("pg14").set("data", "dset6");
    model.result("pg14").create("lngr1", "LineGraph");
    model.result("pg14").feature("lngr1").selection().set(new int[]{1});
    model.result("pg15").set("data", "dset6");
    model.result("pg15").create("lngr1", "LineGraph");
    model.result("pg15").feature("lngr1").selection().set(new int[]{1});
    model.result("pg16").set("data", "dset6");
    model.result("pg16").create("lngr1", "LineGraph");
    model.result("pg16").feature("lngr1").set("xdata", "expr");
    model.result("pg16").feature("lngr1").selection().set(new int[]{1});

    model.component("comp1").probe("pdom1").genResult(null);

    model.study("std2").feature("param").set("pname", new String[]{"grid_factor"});
    model.study("std2").feature("param").set("plistarr", new String[]{"10 50 100"});
    model.study("std2").feature("param").set("punit", new String[]{""});
    model.study("std2").feature("time").set("tlist", "range(0,1,20)");

    model.sol("sol1").attach("std1");
    model.sol("sol1").feature("s1").feature("fc1").set("initstep", 0.01);
    model.sol("sol1").feature("s1").feature("fc1").set("minstep", 1.0E-6);
    model.sol("sol1").feature("s1").feature("fc1").set("maxiter", 50);
    model.sol("sol1").runAll();
    model.sol("sol2").attach("std2");
    model.sol("sol2").feature("v1").set("clist", new String[]{"range(0,1,20)"});
    model.sol("sol2").feature("t1").set("tlist", "range(0,1,20)");
    model.sol("sol2").feature("t1").set("maxorder", 2);
    model.sol("sol2").feature("t1").feature("fc1").set("maxiter", 8);
    model.sol("sol2").feature("t1").feature("fc1").set("damp", 0.9);
    model.sol("sol2").feature("t1").feature("fc1").set("jtech", "once");
    model.sol("sol2").runAll();

    model.batch("p1").set("control", "param");
    model.batch("p1").set("pname", new String[]{"grid_factor"});
    model.batch("p1").set("plistarr", new String[]{"10 50 100"});
    model.batch("p1").set("punit", new String[]{""});
    model.batch("p1").set("err", true);
    model.batch("p1").feature("so1").set("seq", "sol2");
    model.batch("p1").feature("so1").set("psol", "sol5");
    model.batch("p1").feature("so1")
         .set("param", new String[]{"\"grid_factor\",\"10\"", "\"grid_factor\",\"50\"", "\"grid_factor\",\"100\""});
    model.batch("p1").attach("std2");
    model.batch("p1").run();

    model.result().dataset("dset2").label("Probe Solution 2");
    model.result().numerical("pev2").set("looplevelinput", new String[]{"manual"});
    model.result().numerical("pev3").set("unit", new String[]{""});
    model.result("pg8").label("J(t)");
    model.result("pg8").set("solrepresentation", "solnum");
    model.result("pg8").set("xlabel", "Time (s)");
    model.result("pg8").set("ylabel", "Displacement current density, x component (A/m^2), Point: 0");
    model.result("pg8").set("windowtitle", "Probe Plot 1");
    model.result("pg8").set("xlabelactive", false);
    model.result("pg8").set("ylabelactive", false);
    model.result("pg8").feature("tblp1").label("Probe Table Graph 1");
    model.result("pg8").feature("tblp1").set("legend", true);
    model.result("pg9").label("Electric Potential (es)");
    model.result("pg9").feature("lngr1").set("xdataexpr", "x/L_cell");

    return model;
  }

  public static Model run2(Model model) {
    model.result("pg9").feature("lngr1").set("xdataunit", "m");
    model.result("pg9").feature("lngr1").set("xdatadescr", "x/L_cell");
    model.result("pg9").feature("lngr1").set("legend", true);
    model.result("pg9").feature("lngr1").set("resolution", "normal");
    model.result("pg9").feature("lngr2").set("xdataexpr", "x/L_cell");
    model.result("pg9").feature("lngr2").set("xdataunit", "m");
    model.result("pg9").feature("lngr2").set("xdatadescr", "x/L_cell");
    model.result("pg9").feature("lngr2").set("legend", true);
    model.result("pg9").feature("lngr2").set("resolution", "normal");
    model.result("pg9").feature("lngr3").set("expr", "");
    model.result("pg9").feature("lngr3").set("unit", "");
    model.result("pg9").feature("lngr3").set("descr", "");
    model.result("pg9").feature("lngr3").set("resolution", "normal");
    model.result("pg10").label("Concentration (tds)");
    model.result("pg10").set("xlabel", "Arc length");
    model.result("pg10").set("ylabel", "Concentration (mol/m<sup>3</sup>)");
    model.result("pg10").set("xlabelactive", false);
    model.result("pg10").set("ylabelactive", false);
    model.result("pg10").feature("lngr1").label("Line Graph");
    model.result("pg10").feature("lngr1").set("expr", "cp1");
    model.result("pg10").feature("lngr1").set("unit", "mol/m^3");
    model.result("pg10").feature("lngr1").set("descr", "Concentration");
    model.result("pg10").feature("lngr1").set("resolution", "normal");
    model.result("pg11").label("Product Concentrations");
    model.result("pg11").set("looplevelinput", new String[]{"last"});
    model.result("pg11").set("xlabel", "x/L_cell (m)");
    model.result("pg11").set("ylabel", "Concentration (mol/dm<sup>3</sup>)");
    model.result("pg11").set("xlabelactive", false);
    model.result("pg11").set("ylabelactive", false);
    model.result("pg11").feature("lngr1").label("propol");
    model.result("pg11").feature("lngr1").set("expr", "cp1");
    model.result("pg11").feature("lngr1").set("unit", "mol/dm^3");
    model.result("pg11").feature("lngr1").set("descr", "Concentration");
    model.result("pg11").feature("lngr1").set("xdataexpr", "x/L_cell");
    model.result("pg11").feature("lngr1").set("xdataunit", "m");
    model.result("pg11").feature("lngr1").set("xdatadescr", "x/L_cell");
    model.result("pg11").feature("lngr1").set("legend", true);
    model.result("pg11").feature("lngr1").set("legendmethod", "manual");
    model.result("pg11").feature("lngr1").set("legends", new String[]{"propol"});
    model.result("pg11").feature("lngr1").set("resolution", "normal");
    model.result("pg11").feature("lngr2").label("metol");
    model.result("pg11").feature("lngr2").set("expr", "cp3");
    model.result("pg11").feature("lngr2").set("unit", "mol/dm^3");
    model.result("pg11").feature("lngr2").set("descr", "Concentration");
    model.result("pg11").feature("lngr2").set("xdataexpr", "x/L_cell");
    model.result("pg11").feature("lngr2").set("xdataunit", "m");
    model.result("pg11").feature("lngr2").set("xdatadescr", "x/L_cell");
    model.result("pg11").feature("lngr2").set("legend", true);
    model.result("pg11").feature("lngr2").set("legendmethod", "manual");
    model.result("pg11").feature("lngr2").set("legends", new String[]{"metol"});
    model.result("pg11").feature("lngr2").set("resolution", "normal");
    model.result("pg11").feature("lngr3").label("CO");
    model.result("pg11").feature("lngr3").set("expr", "cp5");
    model.result("pg11").feature("lngr3").set("unit", "mol/dm^3");
    model.result("pg11").feature("lngr3").set("descr", "Concentration");
    model.result("pg11").feature("lngr3").set("xdataexpr", "x/L_cell");
    model.result("pg11").feature("lngr3").set("xdataunit", "m");
    model.result("pg11").feature("lngr3").set("xdatadescr", "x/L_cell");
    model.result("pg11").feature("lngr3").set("legend", true);
    model.result("pg11").feature("lngr3").set("legendmethod", "manual");
    model.result("pg11").feature("lngr3").set("legends", new String[]{"CO"});
    model.result("pg11").feature("lngr3").set("resolution", "normal");
    model.result("pg11").feature("lngr4").label("etgly");
    model.result("pg11").feature("lngr4").set("expr", "cp6");
    model.result("pg11").feature("lngr4").set("unit", "mol/dm^3");
    model.result("pg11").feature("lngr4").set("descr", "Concentration");
    model.result("pg11").feature("lngr4").set("xdataexpr", "x/L_cell");
    model.result("pg11").feature("lngr4").set("xdataunit", "m");
    model.result("pg11").feature("lngr4").set("xdatadescr", "x/L_cell");
    model.result("pg11").feature("lngr4").set("legend", true);
    model.result("pg11").feature("lngr4").set("legendmethod", "manual");
    model.result("pg11").feature("lngr4").set("legends", new String[]{"etgly"});
    model.result("pg11").feature("lngr4").set("resolution", "normal");
    model.result("pg11").feature("lngr5").label("allyl");
    model.result("pg11").feature("lngr5").set("expr", "cp7");
    model.result("pg11").feature("lngr5").set("unit", "mol/dm^3");
    model.result("pg11").feature("lngr5").set("descr", "Concentration");
    model.result("pg11").feature("lngr5").set("xdataexpr", "x/L_cell");
    model.result("pg11").feature("lngr5").set("xdataunit", "m");
    model.result("pg11").feature("lngr5").set("xdatadescr", "x/L_cell");
    model.result("pg11").feature("lngr5").set("legend", true);
    model.result("pg11").feature("lngr5").set("legendmethod", "manual");
    model.result("pg11").feature("lngr5").set("legends", new String[]{"allyl"});
    model.result("pg11").feature("lngr5").set("resolution", "normal");
    model.result("pg11").feature("lngr6").label("H2");
    model.result("pg11").feature("lngr6").set("expr", "cp8");
    model.result("pg11").feature("lngr6").set("unit", "mol/dm^3");
    model.result("pg11").feature("lngr6").set("descr", "Concentration");
    model.result("pg11").feature("lngr6").set("xdataexpr", "x/L_cell");
    model.result("pg11").feature("lngr6").set("xdataunit", "m");
    model.result("pg11").feature("lngr6").set("xdatadescr", "x/L_cell");
    model.result("pg11").feature("lngr6").set("legend", true);
    model.result("pg11").feature("lngr6").set("legendmethod", "manual");
    model.result("pg11").feature("lngr6").set("legends", new String[]{"H2"});
    model.result("pg11").feature("lngr6").set("resolution", "normal");
    model.result("pg11").feature("lngr7").label("etol");
    model.result("pg11").feature("lngr7").set("expr", "cp11");
    model.result("pg11").feature("lngr7").set("unit", "mol/dm^3");
    model.result("pg11").feature("lngr7").set("descr", "Concentration");
    model.result("pg11").feature("lngr7").set("xdataexpr", "x/L_cell");
    model.result("pg11").feature("lngr7").set("xdataunit", "m");
    model.result("pg11").feature("lngr7").set("xdatadescr", "x/L_cell");
    model.result("pg11").feature("lngr7").set("legend", true);
    model.result("pg11").feature("lngr7").set("legendmethod", "manual");
    model.result("pg11").feature("lngr7").set("legends", new String[]{"etol"});
    model.result("pg11").feature("lngr7").set("resolution", "normal");
    model.result("pg11").feature("lngr8").label("CH4");
    model.result("pg11").feature("lngr8").set("expr", "cp13");
    model.result("pg11").feature("lngr8").set("unit", "mol/dm^3");
    model.result("pg11").feature("lngr8").set("descr", "Concentration");
    model.result("pg11").feature("lngr8").set("xdataexpr", "x/L_cell");
    model.result("pg11").feature("lngr8").set("xdataunit", "m");
    model.result("pg11").feature("lngr8").set("xdatadescr", "x/L_cell");
    model.result("pg11").feature("lngr8").set("legend", true);
    model.result("pg11").feature("lngr8").set("legendmethod", "manual");
    model.result("pg11").feature("lngr8").set("legends", new String[]{"CH4"});
    model.result("pg11").feature("lngr8").set("resolution", "normal");
    model.result("pg11").feature("lngr9").label("C2H4");
    model.result("pg11").feature("lngr9").set("expr", "cp14");
    model.result("pg11").feature("lngr9").set("unit", "mol/dm^3");
    model.result("pg11").feature("lngr9").set("descr", "Concentration");
    model.result("pg11").feature("lngr9").set("xdataexpr", "x/L_cell");
    model.result("pg11").feature("lngr9").set("xdataunit", "m");
    model.result("pg11").feature("lngr9").set("xdatadescr", "x/L_cell");
    model.result("pg11").feature("lngr9").set("legend", true);
    model.result("pg11").feature("lngr9").set("legendmethod", "manual");
    model.result("pg11").feature("lngr9").set("legends", new String[]{"C2H4"});
    model.result("pg11").feature("lngr9").set("resolution", "normal");
    model.result("pg11").feature("lngr10").label("HCOO-");
    model.result("pg11").feature("lngr10").set("expr", "cp15");
    model.result("pg11").feature("lngr10").set("unit", "mol/dm^3");
    model.result("pg11").feature("lngr10").set("descr", "Concentration");
    model.result("pg11").feature("lngr10").set("xdataexpr", "x/L_cell");
    model.result("pg11").feature("lngr10").set("xdataunit", "m");
    model.result("pg11").feature("lngr10").set("xdatadescr", "x/L_cell");
    model.result("pg11").feature("lngr10").set("legend", true);
    model.result("pg11").feature("lngr10").set("legendmethod", "manual");
    model.result("pg11").feature("lngr10").set("legends", new String[]{"HCOO-"});
    model.result("pg11").feature("lngr10").set("resolution", "normal");
    model.result("pg11").feature("lngr11").label("acet");
    model.result("pg11").feature("lngr11").set("expr", "cp17");
    model.result("pg11").feature("lngr11").set("unit", "mol/dm^3");
    model.result("pg11").feature("lngr11").set("descr", "Concentration");
    model.result("pg11").feature("lngr11").set("xdataexpr", "x/L_cell");
    model.result("pg11").feature("lngr11").set("xdataunit", "m");
    model.result("pg11").feature("lngr11").set("xdatadescr", "x/L_cell");
    model.result("pg11").feature("lngr11").set("legend", true);
    model.result("pg11").feature("lngr11").set("legendmethod", "manual");
    model.result("pg11").feature("lngr11").set("legends", new String[]{"acet"});
    model.result("pg11").feature("lngr11").set("resolution", "normal");
    model.result("pg12").label("Buffer Concentrations");
    model.result("pg12").set("looplevelinput", new String[]{"last", "all"});
    model.result("pg12").set("xlabel", "x/L_cell (m)");
    model.result("pg12").set("ylabel", "Concentration (mol/dm<sup>3</sup>)");
    model.result("pg12").set("ylog", true);
    model.result("pg12").set("xlabelactive", false);
    model.result("pg12").set("ylabelactive", false);
    model.result("pg12").feature("lngr1").label("OH-");
    model.result("pg12").feature("lngr1").set("expr", "cp2");
    model.result("pg12").feature("lngr1").set("unit", "mol/dm^3");
    model.result("pg12").feature("lngr1").set("descr", "Concentration");
    model.result("pg12").feature("lngr1").set("xdataexpr", "x/L_cell");
    model.result("pg12").feature("lngr1").set("xdataunit", "m");
    model.result("pg12").feature("lngr1").set("xdatadescr", "x/L_cell");
    model.result("pg12").feature("lngr1").set("legend", true);
    model.result("pg12").feature("lngr1").set("legendmethod", "manual");
    model.result("pg12").feature("lngr1").set("legends", new String[]{"OH-"});
    model.result("pg12").feature("lngr1").set("resolution", "normal");
    model.result("pg12").feature("lngr2").label("CO2");
    model.result("pg12").feature("lngr2").set("expr", "cp4");
    model.result("pg12").feature("lngr2").set("unit", "mol/dm^3");
    model.result("pg12").feature("lngr2").set("descr", "Concentration");
    model.result("pg12").feature("lngr2").set("xdataexpr", "x/L_cell");
    model.result("pg12").feature("lngr2").set("xdataunit", "m");
    model.result("pg12").feature("lngr2").set("xdatadescr", "x/L_cell");
    model.result("pg12").feature("lngr2").set("legend", true);
    model.result("pg12").feature("lngr2").set("legendmethod", "manual");
    model.result("pg12").feature("lngr2").set("legends", new String[]{"CO2"});
    model.result("pg12").feature("lngr2").set("resolution", "normal");
    model.result("pg12").feature("lngr3").label("K+");
    model.result("pg12").feature("lngr3").set("expr", "cp10");
    model.result("pg12").feature("lngr3").set("unit", "mol/dm^3");
    model.result("pg12").feature("lngr3").set("descr", "Concentration");
    model.result("pg12").feature("lngr3").set("xdataexpr", "x/L_cell");
    model.result("pg12").feature("lngr3").set("xdataunit", "m");
    model.result("pg12").feature("lngr3").set("xdatadescr", "x/L_cell");
    model.result("pg12").feature("lngr3").set("legend", true);
    model.result("pg12").feature("lngr3").set("legendmethod", "manual");
    model.result("pg12").feature("lngr3").set("legends", new String[]{"K+"});
    model.result("pg12").feature("lngr3").set("resolution", "normal");
    model.result("pg12").feature("lngr5").label("HCO3-");
    model.result("pg12").feature("lngr5").set("expr", "cp12");
    model.result("pg12").feature("lngr5").set("unit", "mol/dm^3");
    model.result("pg12").feature("lngr5").set("descr", "Concentration");
    model.result("pg12").feature("lngr5").set("xdataexpr", "x/L_cell");
    model.result("pg12").feature("lngr5").set("xdataunit", "m");
    model.result("pg12").feature("lngr5").set("xdatadescr", "x/L_cell");
    model.result("pg12").feature("lngr5").set("legend", true);
    model.result("pg12").feature("lngr5").set("legendmethod", "manual");
    model.result("pg12").feature("lngr5").set("legends", new String[]{"HCO3-"});
    model.result("pg12").feature("lngr5").set("resolution", "normal");
    model.result("pg12").feature("lngr6").label("CO32-");
    model.result("pg12").feature("lngr6").set("expr", "cp16");
    model.result("pg12").feature("lngr6").set("unit", "mol/dm^3");
    model.result("pg12").feature("lngr6").set("descr", "Concentration");
    model.result("pg12").feature("lngr6").set("xdataexpr", "x/L_cell");
    model.result("pg12").feature("lngr6").set("xdataunit", "m");
    model.result("pg12").feature("lngr6").set("xdatadescr", "x/L_cell");
    model.result("pg12").feature("lngr6").set("legend", true);
    model.result("pg12").feature("lngr6").set("legendmethod", "manual");
    model.result("pg12").feature("lngr6").set("legends", new String[]{"CO32-"});
    model.result("pg12").feature("lngr6").set("resolution", "normal");
    model.result("pg13").label("Fluxes");
    model.result("pg13").set("looplevelinput", new String[]{"last", "all"});
    model.result("pg13").set("xlabel", "x/L_cell (m)");
    model.result("pg13").set("ylabel", "-D2*d(cp2,x)-D2*F_const/RT*cp2*es.Ex (mol/(m<sup>2</sup>*s))");
    model.result("pg13").set("axislimits", true);
    model.result("pg13").set("xmin", 0);
    model.result("pg13").set("xmax", 9.999999747378752E-5);
    model.result("pg13").set("xlabelactive", false);
    model.result("pg13").set("ylabelactive", false);
    model.result("pg13").feature("lngr1").label("Total flux");
    model.result("pg13").feature("lngr1").set("expr", "-D2*d(cp2,x)-D2*F_const/RT*cp2*es.Ex");
    model.result("pg13").feature("lngr1").set("unit", "mol/(m^2*s)");
    model.result("pg13").feature("lngr1").set("descr", "-D2*d(cp2,x)-D2*F_const/RT*cp2*es.Ex");
    model.result("pg13").feature("lngr1").set("xdataexpr", "x/L_cell");
    model.result("pg13").feature("lngr1").set("xdataunit", "m");
    model.result("pg13").feature("lngr1").set("xdatadescr", "x/L_cell");
    model.result("pg13").feature("lngr1").set("linemarker", "cycle");
    model.result("pg13").feature("lngr1").set("markerpos", "datapoints");
    model.result("pg13").feature("lngr1").set("legend", true);
    model.result("pg13").feature("lngr1").set("legendprefix", "total, ");
    model.result("pg13").feature("lngr1").set("resolution", "normal");
    model.result("pg13").feature("lngr2").active(false);
    model.result("pg13").feature("lngr2").set("expr", "D2*cp2*d(cp2,x)");
    model.result("pg13").feature("lngr2").set("unit", "mol^2/(m^5*s)");
    model.result("pg13").feature("lngr2").set("descr", "D2*cp2*d(cp2,x)");
    model.result("pg13").feature("lngr2").set("xdataexpr", "x/L_cell");
    model.result("pg13").feature("lngr2").set("xdataunit", "m");
    model.result("pg13").feature("lngr2").set("xdatadescr", "x/L_cell");
    model.result("pg13").feature("lngr2").set("legend", true);
    model.result("pg13").feature("lngr2").set("legendprefix", "diffusion, ");
    model.result("pg13").feature("lngr2").set("resolution", "normal");
    model.result("pg13").feature("lngr3").active(false);
    model.result("pg13").feature("lngr3").set("expr", "D2*F_const/RT*cp2*es.Ex");
    model.result("pg13").feature("lngr3").set("unit", "mol/(m^2*s)");
    model.result("pg13").feature("lngr3").set("descr", "D2*F_const/RT*cp2*es.Ex");
    model.result("pg13").feature("lngr3").set("xdataexpr", "x/L_cell");
    model.result("pg13").feature("lngr3").set("xdataunit", "m");
    model.result("pg13").feature("lngr3").set("xdatadescr", "x/L_cell");
    model.result("pg13").feature("lngr3").set("legend", true);
    model.result("pg13").feature("lngr3").set("legendprefix", "migration, ");
    model.result("pg13").feature("lngr3").set("resolution", "normal");
    model.result("pg14").label("Electric Potential (es) 1");
    model.result("pg14").set("looplevelinput", new String[]{"last", "all"});
    model.result("pg14").set("xlabel", "Arc length");
    model.result("pg14").set("ylabel", "Electric potential (V)");
    model.result("pg14").set("xlabelactive", false);
    model.result("pg14").set("ylabelactive", false);
    model.result("pg14").feature("lngr1").set("resolution", "normal");
    model.result("pg15").label("Concentration (tds) 1");
    model.result("pg15").set("xlabel", "Arc length");
    model.result("pg15").set("ylabel", "Concentration (mol/m<sup>3</sup>)");
    model.result("pg15").set("xlabelactive", false);
    model.result("pg15").set("ylabelactive", false);
    model.result("pg15").feature("lngr1").label("Line Graph");
    model.result("pg15").feature("lngr1").set("expr", "cp1");
    model.result("pg15").feature("lngr1").set("unit", "mol/m^3");
    model.result("pg15").feature("lngr1").set("descr", "Concentration");
    model.result("pg15").feature("lngr1").set("resolution", "normal");
    model.result("pg16").label("pH");
    model.result("pg16").set("xlabel", "x/L_cell (m)");
    model.result("pg16").set("ylabel", "pH");
    model.result("pg16").set("xlabelactive", false);
    model.result("pg16").set("ylabelactive", false);
    model.result("pg16").feature("lngr1").label("pH");
    model.result("pg16").feature("lngr1").set("expr", "14+log10(cp2/1000.)");
    model.result("pg16").feature("lngr1").set("unit", "");
    model.result("pg16").feature("lngr1").set("descractive", true);
    model.result("pg16").feature("lngr1").set("descr", "pH");
    model.result("pg16").feature("lngr1").set("xdataexpr", "x/L_cell");
    model.result("pg16").feature("lngr1").set("xdataunit", "m");
    model.result("pg16").feature("lngr1").set("xdatadescr", "x/L_cell");
    model.result("pg16").feature("lngr1").set("resolution", "normal");

    return model;
  }

  public static void main(String[] args) {
    Model model = run();
    run2(model);
  }

}
