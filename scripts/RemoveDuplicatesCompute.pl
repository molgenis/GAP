use warnings;
use strict;
use Data::Dumper;
main();

sub main {
	InfoMsg("Commandline : $0 ".join(" ",@ARGV)."\n");
	DuplicateParameterRemoval(@ARGV);	
}

sub DuplicateParameterRemoval{
	for my $shFile (@_){
		InfoMsg("File $shFile does not exist!!!")if(! -e $shFile);
		next if(! -e $shFile);
		my $sh = ShLoader($shFile);
		$sh = ShDupRemove($sh);
		#print Dumper(\$sh);
		
		ShWriter($sh);
	}
}

sub ShLoader {
	my $shFile = shift(@_);
	my $self -> {'file'} = $shFile;
	
	InfoMsg ("Loading sh file ".  $self -> {'file'});
	
	open(my $shHandle, "<", $shFile) or die "cannot open $shFile";
	
	while(<$shHandle>){
		$self -> {'line'} -> {$.} = $_;
		push(@{$self -> {'order'}} , $.);
	}
	close $shHandle;
	
	return $self;
}

sub ShWriter {
	my $self = shift @_;
	
	#open(my $shHandle, ">", $shFile) or die "cannot open $shFile";
	#DebugMsg("Making backups and writing to " . $self -> {'file'});
	
	my $mvCmd = "mv -v " . $self -> {'file'}." " . $self -> {'file'}.".bak_".time;
	
	InfoMsg($mvCmd);
	InfoMsg(`$mvCmd`);
	
	
	open(my $shHandle, ">", $self -> {'file'}) or die "cannot open $self -> {'file'}";
	#InfoMsg(Dumper($self -> {'order'}));
	
	for my $lIndex (@{$self -> {'order'}}){
		print $shHandle $self -> {'line'} -> {$lIndex};
	}
	close $shHandle;
}
sub ShDupRemove {
	my $self = shift @_;
	
	InfoMsg("Removing static array dups ".  $self -> {'file'});
	
	my $ignoreLineLen=300;
	my $fileLen = scalar(keys(%{$self -> {'line'}}));
	
	#DebugMsg("$fileLen");
	
	#copy metadata
	my $selfNew;
	$selfNew -> {'file'} = $self -> {'file'};
	$selfNew -> {'order'} = $self -> {'order'};
	
	#process lines removing weaved array dups in long commands
	for my $lIndex (keys(%{$self -> {'line'}})){
		if(length($self -> {'line'} -> {$lIndex}) < $ignoreLineLen){
			$selfNew -> {'line'} -> {$lIndex} = $self -> {'line'} -> {$lIndex};
		}else{
			$selfNew -> {'line'} -> {$lIndex} = ScanArrayDupsAndRemove($self -> {'line'} -> {$lIndex});
			
		}
	}
	#look for static array declarations a la molgenis 
	$selfNew = ShStaticArrayDeclareRemover($selfNew);
		
	return $selfNew;
	
}

sub ShStaticArrayDeclareRemover {
	my $self = shift @_;
	
	InfoMsg("Removing array declare dups ".  $self -> {'file'});
	
	#copy metadata
	my $selfNew;
	$selfNew -> {'file'} = $self -> {'file'};
	$selfNew -> {'order'} = $self -> {'order'};
	
	for my $lIndex (@{$self -> {'order'}}){
		$selfNew -> {'line'} -> {$lIndex} = $self -> {'line'} -> {$lIndex};
		
		warn $lIndex if not ($selfNew -> {'line'} -> {$lIndex});
		
		if($selfNew -> {'line'} -> {$lIndex} =~ m!^(\w*)\[\d+\]\=(\"[\/\w\._\-\*]*\")$!){
			#this code should collect / merge data from array defention blocks and clean then up
			#if same variable as last > continue collecting
			#if other variable -> cleanup prevous variable by reinjecting previous lines and initiate collect of other variable 
			#$selfNew -> {'line'} -> {$lIndex} =~ m!^(\w*)\[\d+\]\=(\"[\/\w\._\-\*]*\")$!;
			#InfoMsg("call cleanup ") if ($selfNew -> {'declare_last'} ne $1);
			
			if(not(defined($selfNew -> {'declare_last'}))){
				$selfNew -> {'declare_last'} = $1;
				$selfNew -> {'declare'} -> {$1} -> {'variablekey'} = $1;
				$selfNew -> {'declare'} -> {$1} -> {'line_first'} = $lIndex;
				$selfNew -> {'declare'} -> {$1} -> {'line_last'} = $lIndex;
				$selfNew -> {'declare'} -> {$1} -> {'variableval'} -> {$lIndex} = $2;
				$selfNew -> {'declare'} -> {$1} -> {'variableuniq'} -> {$2} ++;
				
			}elsif(defined($selfNew -> {'declare_last'}) && $selfNew -> {'declare_last'} eq $1 && $selfNew -> {'declare'} -> {$1} -> {'line_last'} + 1 == $lIndex){
				$selfNew -> {'declare_last'} = $1;
				$selfNew -> {'declare'} -> {$1} -> {'variablekey'} = $1;
				$selfNew -> {'declare'} -> {$1} -> {'line_last'} = $lIndex;
				#$selfNew -> {'declare'} -> {$1} -> {'line_first'} = $lIndex if (not(defined($selfNew -> {'declare'} -> {'line_first'})));
				$selfNew -> {'declare'} -> {$1} -> {'variableval'} -> {$lIndex} = $2;
				$selfNew -> {'declare'} -> {$1} -> {'variableuniq'} ->  {$2} ++;
			}elsif(defined($selfNew -> {'declare_last'}) && $selfNew -> {'declare_last'} ne $1 || $selfNew -> {'declare'} -> {$1} -> {'line_last'} + 1 != $lIndex){
				#cleanup;
				InfoMsg("Removing array declare dups: call CleanUpDuplicatesArrayDefinition ".  $self -> {'file'});
				
				$selfNew=CleanUpDuplicatesArrayDefinition($selfNew, $selfNew -> {'declare_last'});
				
				
				#reassign vars
				$selfNew -> {'line'} -> {$lIndex} =~ m!^(\w*)\[\d+\]\=(\"[\/\w\._\-\*]*\")$!;
				$selfNew -> {'declare_last'} = $1;
				$selfNew -> {'declare'} -> {$1} -> {'variablekey'} = $1;
				$selfNew -> {'declare'} -> {$1} -> {'line_first'} = $lIndex;
				$selfNew -> {'declare'} -> {$1} -> {'variableval'} -> {$lIndex} = $2;
				$selfNew -> {'declare'} -> {$1} -> {'line_last'} = $lIndex;
				$selfNew -> {'declare'} -> {$1} -> {'variableuniq'} ->  {$2} ++;
			}else{
				die "uncaught error $lIndex";
			}
			
			#$selfNew -> {'declare_last'} = {$1};
			#$selfNew -> {'declare'} -> {$1} -> {'variablekey'} = $1;
			#$selfNew -> {'declare'} -> {$1} -> {'line_last'} = $lIndex;
			#$selfNew -> {'declare'} -> {$1} -> {'line_first'} = $lIndex if (not(defined($selfNew -> {'declare'} -> {'line_first'})));
			#$selfNew -> {'declare'} -> {$1} -> {'variableval'} -> {$lIndex} = $2;
		}
	}
	
	if(defined($selfNew -> {'declare_last'})){
		$selfNew = CleanUpDuplicatesArrayDefinition($selfNew, $selfNew -> {'declare_last'}) 
	}
	
	return $selfNew;
	
}
sub CleanUpDuplicatesArrayDefinition {
	my $self = shift @_;
	
	InfoMsg("Removing duplicates of ".$_[0]);
	
	#for safety count the amount of lines between line first / last
	my $key =shift @_;#arrayname to cleanup
	my $count = 0;
	
	for my $lindex (@{$self -> {'order'}}){
		
		if ($self -> {'declare'} -> {$key} -> {'line_last'} >= $lindex && 
		$self -> {'declare'} -> {$key} -> {'line_first'} <= $lindex){

			$count++;

		}
	}
	
	DebugMsg( "something strange has happend more variablevals then lines to declare them on!" , \$count,scalar(keys(%{$self -> {'declare'} -> {$key} -> {'variableval'}})) )if(scalar(keys(%{$self -> {'declare'} -> {$key} -> {'variableval'}})) > $count);
	
	#make array definitions
	my @array = ();
	my $arrayIndex=0;
	for my $arrayval (keys %{$self -> {'declare'} -> {$key} -> {'variableuniq'}}){
		push(@array, ${key}."\[".$arrayIndex."\]=".$arrayval."\n");			
		$arrayIndex++;
	}
	

	#cleanup by redefing some lines and empty stringing the rest
	for my $lIndex (@{$self -> {'order'}}){
		if($self -> {'declare'} -> {$key} -> {'line_last'} >= $lIndex && $self -> {'declare'} -> {$key} -> {'line_first'} <= $lIndex && scalar(@array) > 0){
			$self -> {'line'} -> {$lIndex} = shift @array;
		}elsif($self -> {'declare'} -> {$key} -> {'line_last'} >= $lIndex && $self -> {'declare'} -> {$key} -> {'line_first'} <= $lIndex && scalar(@array) == 0){
			$self -> {'line'} -> {$lIndex} = "";
		}
	}
	return $self;
	
}
sub ScanArrayDupsAndRemove {
	my $line = shift @_;
	my $lineLen= length($line);
	
	if($line =~ m!\"[\/\w\._\-\*]*\"( \"[\/\w\._\-\*]*\")+!){
		my $arraymatchRes;
		$arraymatchRes -> {'match'} = $&;
		my $index = index($line, $arraymatchRes -> {'match'});
		$arraymatchRes -> {'pre'} = substr($line,0,$index);
		$arraymatchRes -> {'post'} = substr($line,$index+length($arraymatchRes -> {'match'}));
		#InfoMsg(Dumper($arraymatchRes));
		
		my @vals = split(' ',$arraymatchRes -> {'match'});
		my %uniq;
		map{$uniq{$_}++;}(@vals);
		$arraymatchRes -> {'match'} =join(" ", keys (%uniq));
		#InfoMsg(Dumper($arraymatchRes));
				
		$arraymatchRes -> {'post'} = ScanArrayDupsAndRemove( $arraymatchRes -> {'post'});
		
		$line = $arraymatchRes -> {'pre'} . $arraymatchRes -> {'match'} . $arraymatchRes -> {'post'};
	}
	
	#for (my $i = 0; $i < $lineLen; $i++){
	#	my $char = substr($line,$i,1);
	#	
	#}
	
	return $line;
}

sub InfoMsg {
	warn "### ".localtime(time())." ### " . join("\n ### ",@_) . "\n"; 
}
sub DebugMsg {
	warn "### ".localtime(time())." ### " . join("\n ### ",@_) . "\n"; 
	die Dumper (\@_);
}
