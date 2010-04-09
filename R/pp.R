#Prints the percentage done of a loop indexed by i out of stop
pp=function(total,count){
	if(missing(count)){count=evalq(i,envir=parent.frame())}
	if(missing(total)){total=evalq(stop,envir=parent.frame())}
	cat(round(100*(count/total)),"%   \r")
}
