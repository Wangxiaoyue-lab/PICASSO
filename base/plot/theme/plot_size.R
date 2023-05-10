inch2cm <- function(inch) {
    inch*2.54
}

cm2inch <- function(cm) {
    cm/2.54
}

px2cm <- function(px,res){
    inch2cm(px/res)
}