new_default <-
function() {
 NAMESPACEpath<-system.file("NAMESPACE",package="phyreg")
 oldpath<-gsub("NAMESPACE","curparms",x=NAMESPACEpath)
 load(oldpath)
 .defparms<-.curparms
 .defparms$control<-NA
 .defparms$test<-NA
 .defparms$data<-NA
 .defparms$taxmatrix<-NA
 .defparms$phydata<-NA
 .defparms$subset<-NA
 .defparms$heightsdata<-NA
 newpath<-gsub("NAMESPACE","defparms",x=NAMESPACEpath)
 save(.defparms,file=newpath)
      }
