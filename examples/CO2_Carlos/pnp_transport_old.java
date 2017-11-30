import com.comsol.model.*;
import com.comsol.model.util.*;
/*
 *Transport Model as created by Transport Framework
 */
public class pnp_transport {
  public static Model run() {
    Model model = ModelUtil.create("Model");
    model.modelPath("/Users/sringe/software/transport/examples");
    model.label("transport");
    model.comments("");
    model.param().set("D1", "1.3e-09[m^2/s]", "n-propanol Diffusion coefficient");
    model.param().set("D2", "5.273e-09[m^2/s]", "hydroxyl Diffusion coefficient");
    model.param().set("D3", "8.4e-10[m^2/s]", "Methanol Diffusion coefficient");
    model.param().set("D4", "1.91e-09[m^2/s]", "carbon dioxide Diffusion coefficient");
    model.param().set("D5", "2.03e-09[m^2/s]", "carbon monoxide Diffusion coefficient");
    model.param().set("D6", "1.102e-09[m^2/s]", "ethyleneglycol Diffusion coefficient");
    model.param().set("D7", "1.1e-09[m^2/s]", "allylalcohol Diffusion coefficient");
    model.param().set("D8", "4.5e-09[m^2/s]", "hydrogen Diffusion coefficient");
    model.param().set("D9", "0.0[m^2/s]", "unknown Diffusion coefficient");
    model.param().set("D10", "1.957e-09[m^2/s]", "potassium Diffusion coefficient");
    model.param().set("D11", "8.4e-10[m^2/s]", "ethanol Diffusion coefficient");
    model.param().set("D12", "1.185e-09[m^2/s]", "bicarbonate Diffusion coefficient");
    model.param().set("D13", "1.49e-09[m^2/s]", "methane Diffusion coefficient");
    model.param().set("D14", "1.87e-09[m^2/s]", "ethylene Diffusion coefficient");
    model.param().set("D15", "1.454e-09[m^2/s]", "formate Diffusion coefficient");
    model.param().set("D16", "9.23e-10[m^2/s]", "carboxylate Diffusion coefficient");
    model.param().set("D17", "1.089e-09[m^2/s]", "acetate Diffusion coefficient");
    model.param().set("Z1", "0", "n-propanol ionic charge");
    model.param().set("Z2", "-1", "hydroxyl ionic charge");
    model.param().set("Z3", "0", "Methanol ionic charge");
    model.param().set("Z4", "0", "carbon dioxide ionic charge");
    model.param().set("Z5", "0", "carbon monoxide ionic charge");
    model.param().set("Z6", "0", "ethyleneglycol ionic charge");
    model.param().set("Z7", "0", "allylalcohol ionic charge");
    model.param().set("Z8", "0", "hydrogen ionic charge");
    model.param().set("Z9", "0", "unknown ionic charge");
    model.param().set("Z10", "1", "potassium ionic charge");
    model.param().set("Z11", "0", "ethanol ionic charge");
    model.param().set("Z12", "-1", "bicarbonate ionic charge");
    model.param().set("Z13", "0", "methane ionic charge");
    model.param().set("Z14", "0", "ethylene ionic charge");
    model.param().set("Z15", "-1", "formate ionic charge");
    model.param().set("Z16", "-2", "carboxylate ionic charge");
    model.param().set("Z17", "-1", "acetate ionic charge");
    model.param().set("k1", "1.32790476537e-07[mol/m^2/s]", "n-propanol flux");
    model.param().set("k2", "0.00187962930218[mol/m^2/s]", "hydroxyl flux");
    model.param().set("k3", "3.38241305932e-08[mol/m^2/s]", "Methanol flux");
    model.param().set("k4", "-0.000132386711388[mol/m^2/s]", "carbon dioxide flux");
    model.param().set("k5", "2.60588829897e-06[mol/m^2/s]", "carbon monoxide flux");
    model.param().set("k6", "2.39269069282e-08[mol/m^2/s]", "ethyleneglycol flux");
    model.param().set("k7", "0.0[mol/m^2/s]", "allylalcohol flux");
    model.param().set("k8", "0.000464243493372[mol/m^2/s]", "hydrogen flux");
    model.param().set("k9", "0.0[mol/m^2/s]", "unknown flux");
    model.param().set("k10", "0.0[mol/m^2/s]", "potassium flux");
    model.param().set("k11", "2.50704426004e-06[mol/m^2/s]", "ethanol flux");
    model.param().set("k12", "0.0[mol/m^2/s]", "bicarbonate flux");
    model.param().set("k13", "9.21248907348e-05[mol/m^2/s]", "methane flux");
    model.param().set("k14", "1.39220331191e-05[mol/m^2/s]", "ethylene flux");
    model.param().set("k15", "4.16043545663e-06[mol/m^2/s]", "formate flux");
    model.param().set("k16", "0.0[mol/m^2/s]", "carboxylate flux");
    model.param().set("k17", "7.86463827995e-08[mol/m^2/s]", "acetate flux");
    model.param().set("ci1", "0.0[mol/m^3]", "n-propanol bulk concentrations");
    model.param().set("ci2", "6.58543643499e-05[mol/m^3]", "hydroxyl bulk concentrations");
    model.param().set("ci3", "0.0[mol/m^3]", "Methanol bulk concentrations");
    model.param().set("ci4", "34.19[mol/m^3]", "carbon dioxide bulk concentrations");
    model.param().set("ci5", "0.0[mol/m^3]", "carbon monoxide bulk concentrations");
    model.param().set("ci6", "0.0[mol/m^3]", "ethyleneglycol bulk concentrations");
    model.param().set("ci7", "0.0[mol/m^3]", "allylalcohol bulk concentrations");
    model.param().set("ci8", "0.0[mol/m^3]", "hydrogen bulk concentrations");
    model.param().set("ci9", "0.0[mol/m^3]", "unknown bulk concentrations");
    model.param().set("ci10", "100.0[mol/m^3]", "potassium bulk concentrations");
    model.param().set("ci11", "0.0[mol/m^3]", "ethanol bulk concentrations");
    model.param().set("ci12", "99.9692958402[mol/m^3]", "bicarbonate bulk concentrations");
    model.param().set("ci13", "0.0[mol/m^3]", "methane bulk concentrations");
    model.param().set("ci14", "0.0[mol/m^3]", "ethylene bulk concentrations");
    model.param().set("ci15", "0.0[mol/m^3]", "formate bulk concentrations");
    model.param().set("ci16", "0.0307041597553[mol/m^3]", "carboxylate bulk concentrations");
    model.param().set("ci17", "0.0[mol/m^3]", "acetate bulk concentrations");
    model.param().set("L_cell", "7.93e-05", "Cell length");
    model.param().set("epsilon", "lambdaD/L_cell", "Dimensionless Debye length scale");
    model.param().set("T", "298[K]", "Temperature");
    model.param().set("RT", "R_const*T", "Molar gas constant * Temperature");
    model.param().set("delta", "lambdaS/lambdaD", "Dimensionless Stern layer thicknesscM 1[mol/m^3] Metal reference concentration");
    model.param().set("alphac", "0.5", "Cathodic charge transfer coefficient");
    model.param().set("alphaa", "1-alphac", "Anodic charge transfer coefficient");
    model.param().set("V", "0.0 [V]", "Potential");
    model.param().set("lambdaD", "9.60661473344e-10[m]", "Debye length");
    model.param().set("CS", "1.8e-05[F/cm^2]", "Stern layer capacitance");
    model.param().set("eps_r", "78.36", "relative permittivity");
    model.param().set("epsS", "epsilon0_const*2.0", "Stern layer effective permittivity");
    model.param().set("lambdaS", "epsS/CS", "Stern layer thicknessV 0.2[V] Electrode potential");
    model.param().set("k1f", "5.93 [1/s*m^3/mol]", "rate constant: CO2 + OH- -> HCO3-");
    model.param().set("k1r", "0.000133558558559 [1/s]", "rate constant: HCO3- -> CO2 + OH-");
    model.param().set("k2f", "100000.0 [1/s*m^3/mol]", "rate constant: HCO3- + OH- -> CO32- + H2O");
    model.param().set("k2r", "21459.2274678 [1/s*m^3/mol]", "rate constant: CO32- + H2O -> HCO3- + OH-");
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
    model.component("comp1").physics("es").create("fp1", "FloatingPotential", 0);
    model.component("comp1").physics().create("tds", "DilutedSpecies", "geom1");
    model.component("comp1").physics("tds").field("concentration").field("cp1");
    model.component("comp1").physics("tds").field("concentration").component(new String[]{"cp1", "cp2", "cp3", "cp4", "cp5", "cp6", "cp7", "cp8", "cp9", "cp10", "cp11", "cp12", "cp13", "cp14", "cp15", "cp16", "cp17"});
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
    model.component("comp1").physics("tds").feature("cdm1").set("z", new String[][]{{"Z1"}, {"Z2"}, {"Z3"}, {"Z4"}, {"Z5"}, {"Z6"}, {"Z7"}, {"Z8"}, {"Z9"}, {"Z10"}, {"Z11"}, {"Z12"}, {"Z13"}, {"Z14"}, {"Z15"}, {"Z16"}, {"Z17"}});
    model.component("comp1").physics("tds").feature("cdm1").set("minput_temperature", "T");
    model.component("comp1").physics("tds").feature("cdm1").set("D_cp1", new String[][]{{"D1"}, {"0"}, {"0"}, {"0"}, {"D1"}, {"0"}, {"0"}, {"0"}, {"D1"}});
    model.component("comp1").physics("tds").feature("cdm1").set("D_cp2", new String[][]{{"D2"}, {"0"}, {"0"}, {"0"}, {"D2"}, {"0"}, {"0"}, {"0"}, {"D2"}});
    model.component("comp1").physics("tds").feature("cdm1").set("D_cp3", new String[][]{{"D3"}, {"0"}, {"0"}, {"0"}, {"D3"}, {"0"}, {"0"}, {"0"}, {"D3"}});
    model.component("comp1").physics("tds").feature("cdm1").set("D_cp4", new String[][]{{"D4"}, {"0"}, {"0"}, {"0"}, {"D4"}, {"0"}, {"0"}, {"0"}, {"D4"}});
    model.component("comp1").physics("tds").feature("cdm1").set("D_cp5", new String[][]{{"D5"}, {"0"}, {"0"}, {"0"}, {"D5"}, {"0"}, {"0"}, {"0"}, {"D5"}});
    model.component("comp1").physics("tds").feature("cdm1").set("D_cp6", new String[][]{{"D6"}, {"0"}, {"0"}, {"0"}, {"D6"}, {"0"}, {"0"}, {"0"}, {"D6"}});
    model.component("comp1").physics("tds").feature("cdm1").set("D_cp7", new String[][]{{"D7"}, {"0"}, {"0"}, {"0"}, {"D7"}, {"0"}, {"0"}, {"0"}, {"D7"}});
    model.component("comp1").physics("tds").feature("cdm1").set("D_cp8", new String[][]{{"D8"}, {"0"}, {"0"}, {"0"}, {"D8"}, {"0"}, {"0"}, {"0"}, {"D8"}});
    model.component("comp1").physics("tds").feature("cdm1").set("D_cp9", new String[][]{{"D9"}, {"0"}, {"0"}, {"0"}, {"D9"}, {"0"}, {"0"}, {"0"}, {"D9"}});
    model.component("comp1").physics("tds").feature("cdm1").set("D_cp10", new String[][]{{"D10"}, {"0"}, {"0"}, {"0"}, {"D10"}, {"0"}, {"0"}, {"0"}, {"D10"}});
    model.component("comp1").physics("tds").feature("cdm1").set("D_cp11", new String[][]{{"D11"}, {"0"}, {"0"}, {"0"}, {"D11"}, {"0"}, {"0"}, {"0"}, {"D11"}});
    model.component("comp1").physics("tds").feature("cdm1").set("D_cp12", new String[][]{{"D12"}, {"0"}, {"0"}, {"0"}, {"D12"}, {"0"}, {"0"}, {"0"}, {"D12"}});
    model.component("comp1").physics("tds").feature("cdm1").set("D_cp13", new String[][]{{"D13"}, {"0"}, {"0"}, {"0"}, {"D13"}, {"0"}, {"0"}, {"0"}, {"D13"}});
    model.component("comp1").physics("tds").feature("cdm1").set("D_cp14", new String[][]{{"D14"}, {"0"}, {"0"}, {"0"}, {"D14"}, {"0"}, {"0"}, {"0"}, {"D14"}});
    model.component("comp1").physics("tds").feature("cdm1").set("D_cp15", new String[][]{{"D15"}, {"0"}, {"0"}, {"0"}, {"D15"}, {"0"}, {"0"}, {"0"}, {"D15"}});
    model.component("comp1").physics("tds").feature("cdm1").set("D_cp16", new String[][]{{"D16"}, {"0"}, {"0"}, {"0"}, {"D16"}, {"0"}, {"0"}, {"0"}, {"D16"}});
    model.component("comp1").physics("tds").feature("cdm1").set("D_cp17", new String[][]{{"D17"}, {"0"}, {"0"}, {"0"}, {"D17"}, {"0"}, {"0"}, {"0"}, {"D17"}});
/*
 *-> ADD initial concentrations
 */
    model.component("comp1").physics("tds").feature("init1").set("initc", new String[][]{{"ci1"}, {"ci2"}, {"ci3"}, {"ci4"}, {"ci5"}, {"ci6"}, {"ci7"}, {"ci8"}, {"ci9"}, {"ci10"}, {"ci11"}, {"ci12"}, {"ci13"}, {"ci14"}, {"ci15"}, {"ci16"}, {"ci17"}});
/*
 *-> ADD BOUNDARIES
 */
    model.component("comp1").physics("tds").feature("fl1").set("species", new int[][]{{1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}});
/*
 *-> ADD FLUXES
 */
    model.component("comp1").physics("tds").feature("fl1").set("N0", new String[][]{{"k0*flux_factor"}, {"k1*flux_factor"}, {"k2*flux_factor"}, {"k3*flux_factor"}, {"k4*flux_factor"}, {"k5*flux_factor"}, {"k6*flux_factor"}, {"k7*flux_factor"}, {"k8*flux_factor"}, {"k9*flux_factor"}, {"k10*flux_factor"}, {"k11*flux_factor"}, {"k12*flux_factor"}, {"k13*flux_factor"}, {"k14*flux_factor"}, {"k15*flux_factor"}, {"k16*flux_factor"}});
    model.component("comp1").physics("tds").feature("gconstr1").set("constraintExpression", "intop2(cm)-(cref*L)");
    model.component("comp1").physics("tds").feature("gconstr1").active(false);
    model.component("comp1").physics("tds").feature("conc1").set("species", new int[][]{{1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}});
    model.component("comp1").physics("tds").feature("conc1").set("c0", new String[][]{{"ci1"}, {"ci2"}, {"ci3"}, {"ci4"}, {"ci5"}, {"ci6"}, {"ci7"}, {"ci8"}, {"ci9"}, {"ci10"}, {"ci11"}, {"ci12"}, {"ci13"}, {"ci14"}, {"ci15"}, {"ci16"}, {"ci17"}});
/*
 *REACTIONS
 */
    model.component("comp1").physics("tds").feature("reac1").set("R_cp2", "-cp4*cp2*k1f+cp12*k1r-cp12*cp2*k2f+cp16*k2r");
    model.component("comp1").physics("tds").feature("reac1").set("R_cp4", "-cp4*cp2*k1f+cp12*k1r");
    model.component("comp1").physics("tds").feature("reac1").set("R_cp12", "+cp4*cp2*k1f-cp12*k1r-cp12*cp2*k2f+cp16*k2r");
    model.component("comp1").physics("tds").feature("reac1").set("R_cp16", "+cp12*cp2*k2f-cp16*k2r");
    model.component("comp1").physics("tds").feature("reac1").active(false);
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
    model.sol("sol1").feature("s1").create("fc1", "FullyCoupled");
    model.sol("sol1").feature("s1").create("d1", "Direct");
    model.sol("sol1").feature("s1").feature().remove("fcDef");
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
