#ifndef Uintah_Component_Arches_CoalHelper_h
#define Uintah_Component_Arches_CoalHelper_h

#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <vector>
#include <CCA/Components/Arches/ChemMix/ChemHelper.h>

namespace Uintah{

  class CoalHelper {


    public:

      /** @brief Return the CoalHelper **/
      static CoalHelper& self(){
        static CoalHelper s;
        return s;
      }

      /** @brief An object to hold the coal db information **/
      struct CoalDBInfo{
        public:
        int Nenv;
        double value;
        double rhop_o;
        double pi;
        double raw_coal_mf;
        double char_mf;
        double ash_mf;
        double mw_avg;
        double h_c0;
        double h_ch0;
        double h_a0;
        double ksi;
        double Tar_fraction;
        double T_hemisphere;        ///< Ash hemispherical temperature
        double T_fluid;             ///< Ash fluid temperature
        double T_soft;              ///< Ash softening temperature
        double T_porosity;          ///< Ash porosity temperature
        double visc_pre_exponential_factor; ///< Ash viscosity pre-exponential factor [poise/K] -Urbain viscosity model
        double visc_activation_energy; ///< Ash viscosity pre-exponential factor [poise/K] -Urbain viscosity model

        std::vector<double> init_ash;
        std::vector<double> init_rawcoal;
        std::vector<double> init_char;
        std::vector<double> sizes;
        std::vector<double> denom;

        std::string rawcoal_base_name;
        std::string char_base_name;

        struct CoalAnalysis{
          double C;
          double H;
          double O;
          double N;
          double S;
          double CHAR;
          double ASH;
          double H2O;
        };

        CoalAnalysis coal;

      };

      /** @brief Parse coal information for use later **/
      void parse_for_coal_info( ProblemSpecP& db ){

        const ProblemSpecP db_root = db->getRootNode();

        double pi = acos(-1.0);
        if ( db_root->findBlock("CFD")->findBlock("ARCHES")->findBlock("ParticleProperties") ){

          ProblemSpecP db_coal_props = db_root->findBlock("CFD")->findBlock("ARCHES")->findBlock("ParticleProperties");

          std::string particleType;
          db_coal_props->getAttribute("type",particleType);
          if ( particleType == "coal" ){

            db_coal_props->require("density",_coal_db.rhop_o);
            db_coal_props->require("diameter_distribution", _coal_db.sizes);
            db_coal_props->require("raw_coal_enthalpy", _coal_db.h_c0);
            db_coal_props->require("char_enthalpy", _coal_db.h_ch0);
            db_coal_props->require("ash_enthalpy", _coal_db.h_a0);
            db_coal_props->getWithDefault( "ksi",_coal_db.ksi,1.0); // Fraction of the heat released by char oxidation that goes to the particle





            if ( db_coal_props->findBlock("ultimate_analysis")){

              //<!-- as received mass fractions C+H+O+N+S+char+ash+moisture=1 -->
              ProblemSpecP db_ua = db_coal_props->findBlock("ultimate_analysis");

              db_ua->require("C",_coal_db.coal.C);
              db_ua->require("H",_coal_db.coal.H);
              db_ua->require("O",_coal_db.coal.O);
              db_ua->require("N",_coal_db.coal.N);
              db_ua->require("S",_coal_db.coal.S);
              db_ua->require("H2O",_coal_db.coal.H2O);
              db_ua->require("ASH",_coal_db.coal.ASH);
              db_ua->require("CHAR",_coal_db.coal.CHAR);

              double coal_daf = _coal_db.coal.C + _coal_db.coal.H
                + _coal_db.coal.O + _coal_db.coal.N + _coal_db.coal.S; //dry ash free coal
              double coal_dry = _coal_db.coal.C + _coal_db.coal.H
                + _coal_db.coal.O + _coal_db.coal.N + _coal_db.coal.S
                + _coal_db.coal.ASH + _coal_db.coal.CHAR; //moisture free coal
              _coal_db.raw_coal_mf = coal_daf / coal_dry;
              _coal_db.char_mf = _coal_db.coal.CHAR / coal_dry;
              _coal_db.ash_mf = _coal_db.coal.ASH / coal_dry;

              _coal_db.init_char.clear();
              _coal_db.init_rawcoal.clear();
              _coal_db.init_ash.clear();
              _coal_db.denom.clear();

              _coal_db.Nenv = _coal_db.sizes.size();

              for ( unsigned int i = 0; i < _coal_db.sizes.size(); i++ ){

                double mass_dry = (pi/6.0) * pow(_coal_db.sizes[i],3) * _coal_db.rhop_o;     // kg/particle
                _coal_db.init_ash.push_back(mass_dry  * _coal_db.ash_mf);                    // kg_ash/particle (initial)
                _coal_db.init_char.push_back(mass_dry * _coal_db.char_mf);                   // kg_char/particle (initial)
                _coal_db.init_rawcoal.push_back(mass_dry * _coal_db.raw_coal_mf);            // kg_ash/particle (initial)
                _coal_db.denom.push_back( _coal_db.init_ash[i] +
                    _coal_db.init_char[i] +
                    _coal_db.init_rawcoal[i] );

              }
              _coal_db.pi = pi;

              double yElem [5];
              yElem[0]=_coal_db.coal.C/coal_daf; // C daf
              yElem[1]=_coal_db.coal.H/coal_daf; // H daf
              yElem[2]=_coal_db.coal.N/coal_daf; // N daf
              yElem[3]=_coal_db.coal.O/coal_daf; // O daf
              yElem[4]=_coal_db.coal.S/coal_daf; // S daf

              double MW [5] = { 12., 1., 14., 16., 32.}; // Atomic weight of elements (C,H,N,O,S) - kg/kmol
              double mw_avg = 0.0; // Mean atomic weight of coal
              for(int i=0;i<5;i++){
                mw_avg += yElem[i]/MW[i];
              }
              _coal_db.mw_avg = 1.0/mw_avg;

              // ESTIMATEE TAR FRACTION - based on:
              // Alexander Josephson et al. 2018-2019 Reduction of a Detailed Soot
              // Model for Simulation of Pyrolyzing Solid Fuel Combustion
              double gasPressure=101325.; //default
              ChemHelper& helper = ChemHelper::self();
              ChemHelper::TableConstantsMapType the_table_constants = helper.get_table_constants();
              std::cout << the_table_constants->size() << " \n";
              auto press_iter = the_table_constants->find("Pressure");
              if ( press_iter != the_table_constants->end() ){
                gasPressure=press_iter->second;
              }
              const double Oc =_coal_db.coal.O/_coal_db.coal.C*12.011/16.; // Oxygen to Carbon molar ratio
              const double Hc =_coal_db.coal.H/_coal_db.coal.C*12.011/1.008; // Hydrogen to carbon molar ratio 
              double Vol;
              db_coal_props->require("daf_volatiles_fraction",Vol); 
              const double Pres = log10(gasPressure/101325.); // log10( pressure in atmospheres) assume atmospheric
              const double ytar = (-124.2+35.7*Pres+93.5*Oc-223.9*Oc*Oc+284.8*Hc-107.3*Hc*Hc+
                              5.48*Vol+0.014*Vol*Vol-58.2*Pres*Hc-0.521*Pres*Vol-5.32*Hc*Vol)/
                             (-303.8+52.4*Pres+1.55E3*Oc-2.46E3*Oc*Oc+656.9*Hc-266.3*Hc*Hc+15.9*Vol+
                                        0.025*Vol*Vol-90.0*Pres*Hc-462.5*Oc*Hc+4.8*Oc*Vol-17.8*Hc*Vol);

              db_coal_props->getWithDefault( "Tar_fraction",_coal_db.Tar_fraction,ytar);  // user can specify their own Tar_fraction
              proc0cout << " Tar_fraction set to: " <<  _coal_db.Tar_fraction  << " \n";
              if  (ytar < 0.025 || ytar > .975){
                throw ProblemSetupException("Something went terribly wrong with the empirical coal tar model.  Please notify David Lignell immediately at 801-422-1772.", __FILE__, __LINE__);
              }

            } else {
              throw ProblemSetupException("Error: No <ultimate_analysis> found in input file.", __FILE__, __LINE__);
            }
            //  ------------------- end Tar fraction estimation -------------------//

            //Ash temperatures:
            db_coal_props->getWithDefault("ash_hemispherical_temperature", _coal_db.T_hemisphere, -999);
            db_coal_props->getWithDefault("ash_fluid_temperature", _coal_db.T_fluid, -999);
            db_coal_props->getWithDefault("ash_softening_temperature", _coal_db.T_soft, -999);
            db_coal_props->getWithDefault("visc_pre_exponential_factor", _coal_db.visc_pre_exponential_factor, -999);
            db_coal_props->getWithDefault("visc_activation_energy", _coal_db.visc_activation_energy, -999);
            _coal_db.T_porosity = 0.5 * (_coal_db.T_soft + _coal_db.T_fluid);

          }
        }
      }

      CoalDBInfo& get_coal_db(){ return _coal_db; }

    private:

      CoalHelper(){}
      ~CoalHelper(){}

      CoalDBInfo _coal_db;

  };
}
#endif
