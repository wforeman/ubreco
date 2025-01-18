
art_make(	
	  LIB_LIBRARIES
      ${ART_FRAMEWORK_CORE}
			art_Utilities canvas
			${MF_MESSAGELOGGER}
			cetlib cetlib_except
      ${Boost_SYSTEM_LIBRARY}
      MODULE_LIBRARIES 
			larcorealg_Geometry
			larcore_Geometry_Geometry_service
			lardataobj_RecoBase
      lardataobj_MCBase
			lardataobj_AnalysisBase	
			lardata_RecoObjects larreco_Deprecated  
			lardata_Utilities
			larreco_RecoAlg
      larreco_Calorimetry
			larevt_Filters
			nusimdata_SimulationBase
			${LARRECO_LIB}
			${LARDATA_LIB}
      ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			art_Persistency_Common canvas
			art_Persistency_Provenance canvas
			${MF_MESSAGELOGGER}
      ${FHICLCPP}
			cetlib cetlib_except
			${CLHEP}
      ${ROOT_GEOM}
      ${ROOT_XMLIO}
      ${ROOT_GDML}
			${ROOT_BASIC_LIB_LIST}
      ubreco_BlipUtils
      ubreco_BlipRecoAlg

)

add_subdirectory(Utils)
add_subdirectory(Alg)
add_subdirectory(job)
install_fhicl()
install_source()
install_headers()


# install txt files
#file(GLOB txt_file *.txt )
#install( FILES ${txt_file} DESTINATION ${product}/${version}/fcl" )
file(GLOB txt_file *.txt)
install_fw(LIST ${txt_file})
