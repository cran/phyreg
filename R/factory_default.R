factory_default <-
function(){

save_defparms<-function() {
 NAMESPACEpath<-system.file("NAMESPACE",package="phyreg")
 newpath<-gsub("NAMESPACE","defparms",x=NAMESPACEpath)
  save(.defparms,file=newpath) 
 }

.defparms<-list(control=NA, test=NA, subset=NULL,data=NA, phydata=NULL, heightsdata=NULL, addDF=0, rho=-1, lorho=0.3, hirho=0.6, errrho=0.02, minrho=0.0001, tolerance=0.000001, oppf=5, opdf=0, parmx=0, parmxz=0, opfunccall=0, taxmatrix=NULL, linputs=FALSE, sinputs=FALSE, means=FALSE, lmshortx=FALSE, lmshortxz=FALSE, lmlongx=FALSE, lmlongxz=FALSE, hinput=FALSE, paper=FALSE, dfwarning=TRUE, oprho=FALSE, reset=FALSE)

save_defparms()
}
