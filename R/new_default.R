new_default <-
function() {
 defparms<-.accEnv$.curparms
 defparms$control<-NA
 defparms$test<-NA
 defparms$data<-NA
 defparms$taxmatrix<-NA
 defparms$phydata<-NA
 defparms$subset<-NA
 defparms$heightsdata<-NA
 .accEnv$.defparms<-defparms
   }
