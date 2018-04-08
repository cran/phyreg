.__global__ <-
c(".defparms", ".curparms")
.outputstoredinf <-
function(inpr, details, requests, H0model, HAmodel) {

if (missing(details)) if (missing(inpr)) stop("No information supplied to work with") else details<-inpr$details

if (missing(requests)) stop("No requests supplied to fulfil")

if (missing(inpr)) {if (missing(H0model) | missing(HAmodel)) stop (".outputstoredinf needs model details from somewhere!!")} else
  {H0model<-inpr$H0model; HAmodel<-inpr$HAmodel}
    
  if (missing(requests)) {cat(paste('You need to supply requests, e.g. requests=c("means", "opdf"), to this function\n')); stop("Please rerun")}
  
  defaultactionlist<-list(oppf=FALSE, opdf=FALSE,oprho=FALSE, dwarnings=FALSE)

with(details, {

if (!is.null(requests[["oppf"]])) if (requests[["oppf"]]) {
fop<-paste("P = ",format(poff,digits=requests[["oppf"]])," for ","F(",testnumdf,",",testdendf,") = ",format(testf,digits=requests[["oppf"]]),sep="")
fop<-paste(fop," for H0: ",format(H0model),", HA: ",format(HAmodel),sep="")
fop<-paste(fop,"\n",sep="")
cat(fop)}

ndri<-function(x) return(1+trunc(max(sapply(x,log10))))


if(!is.null(requests[["opdf"]])) if(requests[["opdf"]])  {

  iwidth<-1+ndri(c(nspec,nommiss,nspecactive))

  fop<-paste("\n\nBreakdown of species numbers:", "\n\n",sep="")
  fop<-paste(fop,"    Total number of species: ",format(nspec,width=iwidth),"\n","          Omitted by subset: ",format(nomspuse,width=iwidth),"\n",sep="")
  fop<-paste(fop," Omitted for missing values: ", format(nommiss,width=iwidth),"\n",sep="");
  fop<-paste(fop,"       Included in analysis: ",format(nspecactive,width=iwidth),"\n\n\n",sep="") 
  

iwidth<-1+ndri(c(nphyhinodes,hilostomsp,hilostvar,shorttotdf))

  fop<-paste(fop, "Breakdown of numbers: of higher nodes:", "\n\n",sep="")
  fop<-paste(fop,"                                                     Total: ",format(nphyhinodes,width=iwidth),"\n","                              Lost through omitted species: ",format(hilostomsp,width=iwidth),"\n",sep="");
  fop<-paste(fop,"                          Lost through lack of variability: ", format(hilostvar,width=iwidth),"\n",sep="");
  fop<-paste(fop,"Remainder as total degrees of freedom for short regression: ",format(shorttotdf,width=iwidth),sep="") 


iwidth<-1+ndri(c(dflx,dfxlost,dfsx,testlongdf,testlost,testnumdf,shorttotdf,shortcondf,addDF,testnumdf,testdendf))

  fop<-paste(fop, "\n\nBreakdown of degrees of freedom in the short regression:", "\n\n",sep="")
  fop<-paste(fop,"Control variables.    In long regression: ",format(dflx,width=iwidth),"\n",sep="")
  fop<-paste(fop,"           Lost for phylogenetic reasons: ",format(dfxlost,width=iwidth),"\n",sep="")
  fop<-paste(fop,"                    Net, by subtraction:  ",format(dfsx,width=iwidth),"\n\n",sep="")
  fop<-paste(fop,"Test variables.       In long regression: ",format(testlongdf,width=iwidth),"\n",sep="")
  fop<-paste(fop,"           Lost for phylogenetic reasons: ",format(testlost,width=iwidth),"\n",sep="")
  fop<-paste(fop,"                     Net, by subtraction: ",format(testnumdf,width=iwidth),"\n\n",sep="")
  fop<-paste(fop,"Hence, the total degrees of freedom break down as follows:","\n\n",sep="");
  fop<-paste(fop,"                                   Total: ",format(shorttotdf,width=iwidth),"\n",sep="")
  fop<-paste(fop,"                                 Control: ",format(shortcondf,width=iwidth),"\n",sep="")
  fop<-paste(fop," Additional fitted parameters (e.g. rho): ",format(addDF,width=iwidth),"\n",sep="")
  fop<-paste(fop,"                                    Test: ",format(testnumdf,width=iwidth),"\n",sep="")
  fop<-paste(fop,"                Residual, by subtraction: ",format(testdendf,width=iwidth),"\n\n",sep="") 
  fop<-paste(fop,"The last two numbers are the numerator and denominator","\n","degrees of freedom in the F-test of the short regression.\n\n",sep="")
  cat(fop)  }

 if(!is.null(requests[["oprho"]])) if(requests[["oprho"]]) {
   fop<-paste("\n","rho = ",format(rho,digits=requests[["oprho"]]),sep="");
   if(edge==2) fop<-paste(fop,", user-determined\n",sep="") else {
                 fop<-paste(fop," found by maximum likelihood (loglik=", format(lik,digits=requests[["oprho"]]),")\n",sep="")
                 if (edge==1) fop<-paste(fop,"\n","NOTE: likelihood reached its maximum at rho's minimum permitted value: perhaps try reducing minrho?\n",sep="")
           }
  cat(fop)}

if(!is.null(requests[["dfwarning"]])) if (requests[["dfwarning"]]) {
fop<-paste("No degrees of freedom missing for phylogenetic reasons\n");
if (length(missingnodes)>1L || testlost>0L) {
fop<-""
if (length(missingnodes)>1L) {
   onmissingnodes<-on[intersect(missingnodes,(nspec+1):(nspec+nhinodes))]  
   lomn<-length(onmissingnodes)
   fop<-paste(fop, "Note: ", lomn, ' denominator degrees of freedom lost  for phylogenetic reasons',"\n",sep="")
   fop<-paste(fop, '    ID numbers of the missing nodes: ', sep="") # Trying more direct (i.e. not processing onmissingnodes into outomn)
   fop<-paste(fop, onmissingnodes,sep=" ") # Two calls so as to use two different sep's.
   }
if (testlost>0) {if (testlost==1) 
        fop<-paste(fop, '(Note: 1 numerator degree of freedom was lost for phylogenetic reasons)',"\n",sep="") else
        fop<-paste(fop, '(Note: ', testlost,' numerator degrees of freedom were lost for phylogenetic reasons)',"\n",sep="") }
}
cat(fop); cat("\n") }

} ) # end of with(details, and so there is one close parenthesis as well as the brace 

}
