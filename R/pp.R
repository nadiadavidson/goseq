#############################################################################
#Description: Prints progress through a loop
#Notes:
#Author: Matthew Young
#Date Modified: 13/12/2010

pp=function(total,count){
	if(missing(count)){count=evalq(i,envir=parent.frame())}
	if(missing(total)){total=evalq(stop,envir=parent.frame())}
	cat(round(100*(count/total)),"%   \r")
}
