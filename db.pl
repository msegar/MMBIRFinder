#!/usr/bin/perl

use DBI;
use DBD::mysql;

my $start_time = time();

$db = "test";
$host = "localhost";
$database = "test";
$user = "**** ENTER USERNAME HERE ****";
$pw = "**** ENTER PASSWORD HERE ****";

# assign the values to your connection variable
$connectionInfo="dbi:mysql:$db;$host";

# make connection to database
$connection = DBI->connect($connectionInfo,$user,$pw);

#READ IN UNALIGNED FILE
#my $filename = "cons_locs.txt";
my $filename = "half_read_clusters.txt";
open(my $fh, $filename) or die "Could not open half_read_clusters.txt file $!";
open (MYFILE, '>mysql_results.txt') or die "Could not open mysql_results.txt file $!";

$counter = 0;
$misses = 0;
$hits = 0;

while (my $line = <$fh>) {
	if ($counter % 10000 == 0) {
		print $counter."\n";
		#print $counter." of 1397028\n";
	}
	$counter += 1;
	
	#Parse the data
	chomp $line;
	
	my @line = split("\t", $line);
	
	#$read = "`".$line[0]."`";
	$read = $line[0];
	$seq = $line[1];
	$loc = $line[2];
	$chr = $line[3];
	$cluster = $line[4];
	$read_orig = $read;
	
	#print "Read orig: ".$read."\n";
	
	$last = substr $read,-1,1;
	
	if ($last == 1) {
		substr ($read, -1, 1, "2");
	} 
	elsif ($last == 2) {
		substr ($read, -1, 1, "1");
	} 
	
	#$input = <STDIN>;
	#print "Read after: ".$read."\n";
	
	$query = "SELECT * FROM `unaligned` WHERE `read`='$read'";
	#print "query = ".$query."\n";
	$statement = $connection->prepare($query);
	$statement->execute();
	#print "size: ".$statement->rows()."\n";
	
	if (@data = $statement->fetchrow_array()){
		$hits = $hits + 1;
		$seq2 = @data[1];
		
		if ($last == 2) {
			$len = length($seq2);
			#print "End 2: \n";
			$loc = $loc - length($seq2);
			print MYFILE $seq2.$seq."\t";
			print MYFILE $read_orig."\t";
			print MYFILE @data[0]."\t";
			#print $read_orig."\t";
			#print @data[0]."\t";
			print MYFILE ($loc)."\t";
			#print ($loc)."\t";
			print MYFILE $chr."\t";
			print MYFILE $cluster."\n";
		}
		elsif ($last == 1) {
			#print "End 1: \n";
			print MYFILE $seq.$seq2."\t";
			print MYFILE $read_orig."\t";
			print MYFILE @data[0]."\t";
			#print $read_orig."\t";
			#print @data[0]."\t";
			print MYFILE $loc."\t";
			#print $loc."\t";
			print MYFILE $chr."\t";
			print MYFILE $cluster."\n";
		}
		else {
			print "uh oh. There was an error\n";
			print "read: ".$read."\n";
			print "last: ".$last."\n";
			print "seq: ".$seq."\n";
			print "seq2: ".$seq."\n";
		}
	} else {
		#print "miss\n";
		$misses = $misses + 1;
	}
	
	#$input = <STDIN>;
	#chomp($input);
}
close(my $filename);
close(MYFILE);

$statement->finish();
$connection->disconnect();

# Print the final time
my $end_time = time();
my $run_time = $end_time - $start_time;
print "Job took $run_time seconds\n";
print "Hits: ".$hits."\n";
print "Misses: ".$misses."\n";
