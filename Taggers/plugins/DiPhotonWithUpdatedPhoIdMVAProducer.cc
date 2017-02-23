#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "flashgg/MicroAOD/interface/PhotonIdUtils.h"
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"

#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandFlat.h"

#include "TFile.h"
#include "TGraph.h"
#include "TVectorD.h"

using namespace std;
using namespace edm;
using namespace reco;

namespace flashgg {

    class DiPhotonWithUpdatedPhoIdMVAProducer : public edm::EDProducer
    {
    public:
        DiPhotonWithUpdatedPhoIdMVAProducer( const edm::ParameterSet & );
        void produce( edm::Event &, const edm::EventSetup & ) override;

    private:
        /** random number generator for simple photon iso smearing */
        double simpleSmearPhotonIso(edm::Event &event, double eta, double iso);

        edm::EDGetTokenT<edm::View<flashgg::DiPhotonCandidate> > token_;
        edm::EDGetTokenT<double> rhoToken_;
        PhotonIdUtils phoTools_;
        edm::FileInPath phoIdMVAweightfileEB_, phoIdMVAweightfileEE_, correctionFile_, non5x5correctionFile_;
        bool correctInputs_;
        bool debug_;
        //        std::vector<TGraph*> corrections_;
        std::vector<std::unique_ptr<TGraph> > corrections_;

        bool doNon5x5transformation_;
        std::vector<std::unique_ptr<TGraph> > non5x5corrections_;

        bool useNewPhoId_;

        EffectiveAreas _effectiveAreas;
        vector<double> _phoIsoPtScalingCoeff;
        double _phoIsoCutoff;
        bool _doSimpleIsoCorrection;

        // index 0 is barrel, 1 is endcap
        TGraph *simplePhoIsoTransformation[2];

        // fraction of events at photon isolation zero
        // index 0 is barrel, 1 is endcap
        double simplePhoIsoZeroFraction_MC[2], simplePhoIsoZeroFraction_data[2]; 
        
        edm::Service<edm::RandomNumberGenerator> simplePhotonIsoRNG;
    };

    DiPhotonWithUpdatedPhoIdMVAProducer::DiPhotonWithUpdatedPhoIdMVAProducer( const edm::ParameterSet &ps ) :
        token_(consumes<edm::View<flashgg::DiPhotonCandidate> >(ps.getParameter<edm::InputTag>("src"))),
        rhoToken_( consumes<double>( ps.getParameter<edm::InputTag>( "rhoFixedGridCollection" ) ) ),
        debug_( ps.getParameter<bool>( "Debug" ) ),
        _effectiveAreas((ps.getParameter<edm::FileInPath>("effAreasConfigFile")).fullPath()),
        _phoIsoPtScalingCoeff(ps.getParameter<std::vector<double >>("phoIsoPtScalingCoeff")),
        _phoIsoCutoff(ps.getParameter<double>("phoIsoCutoff")),
        _doSimpleIsoCorrection(ps.getParameter<bool>("doSimpleIsoCorrection"))
    {
        if (_doSimpleIsoCorrection) {
            cout << "applying simple photon isolation smearing" << endl;
            string fname = ps.getParameter<edm::FileInPath>("simpleIsoCorrectionFile").fullPath();

            // get the correction graphs for barrel and endcap
            // make a clone of the graph so we can close the file
            TFile *fin = new TFile(fname.c_str());
            
            // get the fractions of data and MC at zero bin
            simplePhoIsoTransformation[0] = (TGraph*)(fin->Get("transfPhIsoEB")->Clone());
            simplePhoIsoTransformation[1] = (TGraph*)(fin->Get("transfPhIsoEE")->Clone());

            TVectorD *vec;
            vec = (TVectorD*) fin->Get("phoIsoZeroFraction_MC_EB");
            simplePhoIsoZeroFraction_MC[0]   = (*vec)[0];
            vec = (TVectorD*) fin->Get("phoIsoZeroFraction_data_EB");
            simplePhoIsoZeroFraction_data[0] = (*vec)[0];
            vec = (TVectorD*) fin->Get("phoIsoZeroFraction_MC_EE");
            simplePhoIsoZeroFraction_MC[1]   = (*vec)[0];
            vec = (TVectorD*) fin->Get("phoIsoZeroFraction_data_EE");
            simplePhoIsoZeroFraction_data[1] = (*vec)[0];

            fin->Close();
            delete fin;

            if( ! simplePhotonIsoRNG.isAvailable() ) {
                throw cms::Exception( "Configuration" ) << "simple photon isolation smearing requires the RandomNumberGeneratorService  - please add to configuration";
            }

        } else {
            for (unsigned index = 0; index < 2; ++index) {
                simplePhoIsoTransformation[index]    = NULL;
                simplePhoIsoZeroFraction_MC[index]   = -1;
                simplePhoIsoZeroFraction_data[index] = -1;
            }
        }

        useNewPhoId_ = ps.getParameter<bool>( "useNewPhoId" );
        phoIdMVAweightfileEB_ = ps.getParameter<edm::FileInPath>( "photonIdMVAweightfile_EB" );
        phoIdMVAweightfileEE_ = ps.getParameter<edm::FileInPath>( "photonIdMVAweightfile_EE" );
        if(useNewPhoId_){
            phoIdMVAweightfileEB_ = ps.getParameter<edm::FileInPath>( "photonIdMVAweightfile_EB_new" );
            phoIdMVAweightfileEE_ = ps.getParameter<edm::FileInPath>( "photonIdMVAweightfile_EE_new" );
        }
        phoTools_.setupMVA( phoIdMVAweightfileEB_.fullPath(), phoIdMVAweightfileEE_.fullPath(), useNewPhoId_ );

        correctInputs_ = ps.existsAs<edm::FileInPath>("correctionFile") ? true: false;
        if (correctInputs_) {
            correctionFile_ = ps.getParameter<edm::FileInPath>( "correctionFile" );
            TFile* f = TFile::Open(correctionFile_.fullPath().c_str());
            corrections_.emplace_back((TGraph*)((TGraph*) f->Get("transffull5x5R9EB"))->Clone() );
            corrections_.emplace_back((TGraph*)((TGraph*) f->Get("transfEtaWidthEB"))->Clone() );
            corrections_.emplace_back((TGraph*)((TGraph*) f->Get("transfS4EB"))->Clone() );
            corrections_.emplace_back((TGraph*)((TGraph*) f->Get("transffull5x5R9EE"))->Clone() );
            corrections_.emplace_back((TGraph*)((TGraph*) f->Get("transfEtaWidthEE"))->Clone() );
            corrections_.emplace_back((TGraph*)((TGraph*) f->Get("transfS4EE"))->Clone() );
            corrections_.emplace_back((TGraph*)((TGraph*) f->Get("transffull5x5sieieEB"))->Clone() );
            corrections_.emplace_back((TGraph*)((TGraph*) f->Get("transffull5x5sieieEE"))->Clone() );
            f->Close();
        }

        doNon5x5transformation_ =ps.getParameter<bool>( "doNon5x5transformation" );
        if (doNon5x5transformation_) {
            non5x5correctionFile_ = ps.getParameter<edm::FileInPath>( "non5x5correctionFile" );
            TFile* non5x5_f = TFile::Open(non5x5correctionFile_.fullPath().c_str());
            non5x5corrections_.emplace_back((TGraph*)((TGraph*) non5x5_f->Get("transfr9EB"))->Clone() );
            non5x5corrections_.emplace_back((TGraph*)((TGraph*) non5x5_f->Get("transfr9EE"))->Clone() );
            non5x5corrections_.emplace_back((TGraph*)((TGraph*) non5x5_f->Get("transfsieieEB"))->Clone() );
            non5x5corrections_.emplace_back((TGraph*)((TGraph*) non5x5_f->Get("transfsieieEE"))->Clone() );
            non5x5corrections_.emplace_back((TGraph*)((TGraph*) non5x5_f->Get("transfsipipEB"))->Clone() );
            non5x5corrections_.emplace_back((TGraph*)((TGraph*) non5x5_f->Get("transfsipipEE"))->Clone() );
            non5x5corrections_.emplace_back((TGraph*)((TGraph*) non5x5_f->Get("transfsieipEB"))->Clone() );
            non5x5corrections_.emplace_back((TGraph*)((TGraph*) non5x5_f->Get("transfsieipEE"))->Clone() );
            non5x5_f->Close();
        }

        produces<std::vector<flashgg::DiPhotonCandidate> >();
    }

    void DiPhotonWithUpdatedPhoIdMVAProducer::produce( edm::Event &evt, const edm::EventSetup & )
    {
        edm::Handle<edm::View<flashgg::DiPhotonCandidate> > objects;
        evt.getByToken( token_, objects );

        edm::Handle<double> rhoHandle;
        evt.getByToken( rhoToken_, rhoHandle );
        const double rhoFixedGrd = *( rhoHandle.product() );

        auto_ptr<std::vector<flashgg::DiPhotonCandidate> > out_obj( new std::vector<flashgg::DiPhotonCandidate>() );

        for (const auto & obj : *objects) {
            flashgg::DiPhotonCandidate *new_obj = obj.clone();
            new_obj->makePhotonsPersistent();
            double leadCorrectedEtaWidth = 0., subLeadCorrectedEtaWidth = 0.;
            if (not evt.isRealData() and correctInputs_) { 
                if (new_obj->getLeadingPhoton().isEB()) {
                    if (this->debug_) {
                        std::cout << new_obj->getLeadingPhoton().full5x5_r9() << std::endl;
                        std::cout << new_obj->getLeadingPhoton().old_r9() << std::endl;
                    }
                    reco::Photon::ShowerShape newShowerShapes = new_obj->getLeadingPhoton().full5x5_showerShapeVariables();
                    newShowerShapes.e3x3 = corrections_[0]->Eval(new_obj->getLeadingPhoton().full5x5_r9())*new_obj->getLeadingPhoton().superCluster()->rawEnergy();
                    newShowerShapes.sigmaIetaIeta = corrections_[6]->Eval(new_obj->getLeadingPhoton().full5x5_sigmaIetaIeta());
                    new_obj->getLeadingPhoton().full5x5_setShowerShapeVariables(newShowerShapes);
                    leadCorrectedEtaWidth = corrections_[1]->Eval(new_obj->getLeadingPhoton().superCluster()->etaWidth());
                    new_obj->getLeadingPhoton().getSuperCluster()->setEtaWidth(leadCorrectedEtaWidth);
                    new_obj->getLeadingPhoton().setS4(corrections_[2]->Eval(new_obj->getLeadingPhoton().s4()));

                    if (this->debug_) {
                        std::cout << new_obj->getLeadingPhoton().full5x5_r9() << std::endl;
                        std::cout << new_obj->getLeadingPhoton().old_r9() << std::endl;
                    }
                }
                
                if (new_obj->getSubLeadingPhoton().isEB()) {
                    reco::Photon::ShowerShape newShowerShapes = new_obj->getSubLeadingPhoton().full5x5_showerShapeVariables();
                    newShowerShapes.e3x3 = corrections_[0]->Eval(new_obj->getSubLeadingPhoton().full5x5_r9())*new_obj->getSubLeadingPhoton().superCluster()->rawEnergy();
                    newShowerShapes.sigmaIetaIeta = corrections_[6]->Eval(new_obj->getSubLeadingPhoton().full5x5_sigmaIetaIeta());
                    new_obj->getSubLeadingPhoton().full5x5_setShowerShapeVariables(newShowerShapes);
                    subLeadCorrectedEtaWidth = corrections_[1]->Eval(new_obj->getSubLeadingPhoton().superCluster()->etaWidth());
                    new_obj->getSubLeadingPhoton().getSuperCluster()->setEtaWidth(subLeadCorrectedEtaWidth);
                    new_obj->getSubLeadingPhoton().setS4(corrections_[2]->Eval(new_obj->getSubLeadingPhoton().s4()));
                }

                if (new_obj->getLeadingPhoton().isEE()) {
                    if (this->debug_) {
                        std::cout << new_obj->getLeadingPhoton().full5x5_r9() << std::endl;
                        std::cout << new_obj->getLeadingPhoton().old_r9() << std::endl;
                    }
                    reco::Photon::ShowerShape newShowerShapes = new_obj->getLeadingPhoton().full5x5_showerShapeVariables();
                    newShowerShapes.e3x3 = corrections_[3]->Eval(new_obj->getLeadingPhoton().full5x5_r9())*new_obj->getLeadingPhoton().superCluster()->rawEnergy();
                    newShowerShapes.sigmaIetaIeta = corrections_[7]->Eval(new_obj->getLeadingPhoton().full5x5_sigmaIetaIeta());
                    new_obj->getLeadingPhoton().full5x5_setShowerShapeVariables(newShowerShapes);
                    leadCorrectedEtaWidth = corrections_[4]->Eval(new_obj->getLeadingPhoton().superCluster()->etaWidth());
                    new_obj->getLeadingPhoton().getSuperCluster()->setEtaWidth(leadCorrectedEtaWidth);
                    new_obj->getLeadingPhoton().setS4(corrections_[5]->Eval(new_obj->getLeadingPhoton().s4()));

                    if (this->debug_) {
                        std::cout << new_obj->getLeadingPhoton().full5x5_r9() << std::endl;
                        std::cout << new_obj->getLeadingPhoton().old_r9() << std::endl;
                    }
                }
                
                if (new_obj->getSubLeadingPhoton().isEE()) {
                    reco::Photon::ShowerShape newShowerShapes = new_obj->getSubLeadingPhoton().full5x5_showerShapeVariables();
                    newShowerShapes.e3x3 = corrections_[3]->Eval(new_obj->getSubLeadingPhoton().full5x5_r9())*new_obj->getSubLeadingPhoton().superCluster()->rawEnergy();
                    newShowerShapes.sigmaIetaIeta = corrections_[7]->Eval(new_obj->getSubLeadingPhoton().full5x5_sigmaIetaIeta());
                    new_obj->getSubLeadingPhoton().full5x5_setShowerShapeVariables(newShowerShapes);
                    subLeadCorrectedEtaWidth = corrections_[4]->Eval(new_obj->getSubLeadingPhoton().superCluster()->etaWidth());
                    new_obj->getSubLeadingPhoton().getSuperCluster()->setEtaWidth(subLeadCorrectedEtaWidth);
                    new_obj->getSubLeadingPhoton().setS4(corrections_[5]->Eval(new_obj->getSubLeadingPhoton().s4()));
                }
            }

            if (not evt.isRealData() and doNon5x5transformation_) { 
                if (new_obj->getLeadingPhoton().isEB()) {
                    if (this->debug_) {
                        std::cout << new_obj->getLeadingPhoton().full5x5_r9() << std::endl;
                        std::cout << new_obj->getLeadingPhoton().old_r9() << std::endl;
                    }
                    reco::Photon::ShowerShape non5x5ShowerShapes = new_obj->getLeadingPhoton().showerShapeVariables();
                    non5x5ShowerShapes.e3x3 = non5x5corrections_[0]->Eval(new_obj->getLeadingPhoton().old_r9())*new_obj->getLeadingPhoton().superCluster()->rawEnergy();
                    non5x5ShowerShapes.sigmaIetaIeta = non5x5corrections_[2]->Eval(new_obj->getLeadingPhoton().sigmaIetaIeta());
                    new_obj->getLeadingPhoton().setShowerShapeVariables(non5x5ShowerShapes);
                    new_obj->getLeadingPhoton().setSipip(non5x5corrections_[4]->Eval(new_obj->getLeadingPhoton().sipip()));
                    new_obj->getLeadingPhoton().setSieip(non5x5corrections_[6]->Eval(new_obj->getLeadingPhoton().sieip()));

                    if (this->debug_) {
                        std::cout << new_obj->getLeadingPhoton().full5x5_r9() << std::endl;
                        std::cout << new_obj->getLeadingPhoton().old_r9() << std::endl;
                    }
                }
                
                if (new_obj->getSubLeadingPhoton().isEB()) {
                    reco::Photon::ShowerShape non5x5ShowerShapes = new_obj->getSubLeadingPhoton().showerShapeVariables();
                    non5x5ShowerShapes.e3x3 = non5x5corrections_[0]->Eval(new_obj->getSubLeadingPhoton().old_r9())*new_obj->getSubLeadingPhoton().superCluster()->rawEnergy();
                    non5x5ShowerShapes.sigmaIetaIeta = non5x5corrections_[2]->Eval(new_obj->getSubLeadingPhoton().sigmaIetaIeta());
                    new_obj->getSubLeadingPhoton().setShowerShapeVariables(non5x5ShowerShapes);
                    new_obj->getSubLeadingPhoton().setSipip(non5x5corrections_[4]->Eval(new_obj->getSubLeadingPhoton().sipip()));
                    new_obj->getSubLeadingPhoton().setSieip(non5x5corrections_[6]->Eval(new_obj->getSubLeadingPhoton().sieip()));
                }

                if (new_obj->getLeadingPhoton().isEE()) {
                    if (this->debug_) {
                        std::cout << new_obj->getLeadingPhoton().full5x5_r9() << std::endl;
                        std::cout << new_obj->getLeadingPhoton().old_r9() << std::endl;
                    }
                    reco::Photon::ShowerShape non5x5ShowerShapes = new_obj->getLeadingPhoton().showerShapeVariables();
                    non5x5ShowerShapes.e3x3 = non5x5corrections_[1]->Eval(new_obj->getLeadingPhoton().old_r9())*new_obj->getLeadingPhoton().superCluster()->rawEnergy();
                    non5x5ShowerShapes.sigmaIetaIeta = non5x5corrections_[3]->Eval(new_obj->getLeadingPhoton().sigmaIetaIeta());
                    new_obj->getLeadingPhoton().setShowerShapeVariables(non5x5ShowerShapes);
                    new_obj->getLeadingPhoton().setSipip(non5x5corrections_[5]->Eval(new_obj->getLeadingPhoton().sipip()));
                    new_obj->getLeadingPhoton().setSieip(non5x5corrections_[7]->Eval(new_obj->getLeadingPhoton().sieip()));

                    if (this->debug_) {
                        std::cout << new_obj->getLeadingPhoton().full5x5_r9() << std::endl;
                        std::cout << new_obj->getLeadingPhoton().old_r9() << std::endl;
                    }
                }
                
                if (new_obj->getSubLeadingPhoton().isEE()) {
                    reco::Photon::ShowerShape non5x5ShowerShapes = new_obj->getSubLeadingPhoton().showerShapeVariables();
                    non5x5ShowerShapes.e3x3 = non5x5corrections_[1]->Eval(new_obj->getSubLeadingPhoton().old_r9())*new_obj->getSubLeadingPhoton().superCluster()->rawEnergy();
                    non5x5ShowerShapes.sigmaIetaIeta = non5x5corrections_[3]->Eval(new_obj->getSubLeadingPhoton().sigmaIetaIeta());
                    new_obj->getSubLeadingPhoton().setShowerShapeVariables(non5x5ShowerShapes);
                    new_obj->getSubLeadingPhoton().setSipip(non5x5corrections_[5]->Eval(new_obj->getSubLeadingPhoton().sipip()));
                    new_obj->getSubLeadingPhoton().setSieip(non5x5corrections_[7]->Eval(new_obj->getSubLeadingPhoton().sieip()));
                }
            }

            //----------
            if (not evt.isRealData() and  _doSimpleIsoCorrection ) {
                float lead_iso = new_obj->getLeadingPhoton().pfPhoIso03();
                float sublead_iso = new_obj->getSubLeadingPhoton().pfPhoIso03();
                float lead_eta = new_obj->getLeadingPhoton().superCluster()->eta();
                float sublead_eta = new_obj->getSubLeadingPhoton().superCluster()->eta();
                if (this->debug_) {
                    std::cout << "Doing Iso correction to lead (sublead) photon with eta,iso: " << lead_eta << ", " << lead_iso;
                    std::cout << " (" << sublead_eta << ", " << ", " << sublead_iso << ")" << std::endl;
                }

                float new_iso_lead = simpleSmearPhotonIso(evt, lead_eta, lead_iso);
                float new_iso_sublead = simpleSmearPhotonIso(evt, sublead_eta, sublead_iso);
                new_obj->getLeadingPhoton().setpfPhoIso03(new_iso_lead);
                new_obj->getSubLeadingPhoton().setpfPhoIso03(new_iso_sublead);
                if (this->debug_) {
                    std::cout << " Final iso value for lead (sublead) photon: " << new_obj->getLeadingPhoton().pfPhoIso03() << " (" 
                              << new_obj->getSubLeadingPhoton().pfPhoIso03() << ")" << std::endl;
                }
            }

            //----------

            if (this->debug_) {
                std::cout << " Input DiPhoton lead (sublead) MVA: " << obj.leadPhotonId() << " " << obj.subLeadPhotonId() << std::endl;
            }
            double eA_leadPho = _effectiveAreas.getEffectiveArea( abs(new_obj->getLeadingPhoton().superCluster()->eta()) );
            double eA_subLeadPho = _effectiveAreas.getEffectiveArea( abs(new_obj->getSubLeadingPhoton().superCluster()->eta()) );
            
            float newleadmva = phoTools_.computeMVAWrtVtx( new_obj->getLeadingPhoton(), new_obj->vtx(), rhoFixedGrd, leadCorrectedEtaWidth, eA_leadPho, _phoIsoPtScalingCoeff, _phoIsoCutoff );
            new_obj->getLeadingPhoton().setPhoIdMvaWrtVtx( new_obj->vtx(), newleadmva);
            float newsubleadmva = phoTools_.computeMVAWrtVtx( new_obj->getSubLeadingPhoton(), new_obj->vtx(), rhoFixedGrd, subLeadCorrectedEtaWidth,eA_subLeadPho, _phoIsoPtScalingCoeff, _phoIsoCutoff );
            new_obj->getSubLeadingPhoton().setPhoIdMvaWrtVtx( new_obj->vtx(), newsubleadmva);
            if (this->debug_) {
                std::cout << " Output DiPhoton lead (sublead) MVA: " << new_obj->leadPhotonId() << " " << new_obj->subLeadPhotonId() << std::endl;
            }
            out_obj->push_back(*new_obj);
            delete new_obj;
        }
        evt.put(out_obj);
    }

    //----------------------------------------------------------------------

    /** do simple photon iso smearing on MC (this should not be called on data) */
    double DiPhotonWithUpdatedPhoIdMVAProducer::simpleSmearPhotonIso(edm::Event &event, double eta, double iso) {

        double abseta = fabs(eta);

        unsigned index;

        if (abseta < 1.5) {
            // barrel
            index = 0;
        } else {
            // endcap
            index = 1;
        }
        
        // no data beween zero and this value
        const double phIsoMinVal = 0.03;

        // range to which to smear the excess MC values at zero
        const double smearBegin = 0.01;
        const double smearEnd = phIsoMinVal;

        CLHEP::HepRandomEngine & engine = simplePhotonIsoRNG->getEngine( event.streamID() );


        if (iso <= phIsoMinVal) {
            // decide whether this an excess event or not
            if (CLHEP::RandFlat::shoot(&engine, 0., 1.) * simplePhoIsoZeroFraction_MC[index] > simplePhoIsoZeroFraction_data[index]) {

                // spread the event evenly to the target range
                iso = CLHEP::RandFlat::shoot(&engine, smearBegin, smearEnd);

                // now apply the transformation
                iso = simplePhoIsoTransformation[index]->Eval(iso);

            } else {
                // event stays at zero, do NOT feed through the correction graph
            }
        } else {
            // non-zero value, apply the correction graph
            iso = simplePhoIsoTransformation[index]->Eval(iso);
        }

        return iso;
    }

    //----------------------------------------------------------------------
}

typedef flashgg::DiPhotonWithUpdatedPhoIdMVAProducer FlashggDiPhotonWithUpdatedPhoIdMVAProducer;
DEFINE_FWK_MODULE( FlashggDiPhotonWithUpdatedPhoIdMVAProducer );

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
