#! /usr/bin/ruby -w

class Array

	def all_equal_to( item )
		self.select{ |e| e == item }.size == self.size
	end

end

class Site

	@@num = 0

	attr_reader :history
	attr_reader :status

	def initialize( status )
		@@num += 1
		@num = @@num
		@current_status = status
		@ref_status = status
		@history = Array.new( @@memory_length, status )
	end

	def update( new_status )
		unless( @@memory_length == 0 )	
			@history << @current_status
			@history.shift
		end
		@current_status = new_status
		@ref_status = @current_status if is_unchanged? or @ref_status.nil?
	end

	def is_unchanged?
		@history.all_equal_to( @current_status )
	end

	def is_a_vacancy?
		@ref_status == 0
	end

	def self.memory_length=( persistence_length )
		@@memory_length = persistence_length
	end

	def remove_vibrations # if a site has returned to its reference state within the persistence time, we can ignore the intermediate steps.
		if ( @ref_status == @current_status ) 
			@history = Array.new(@history.size, @status)
		end
	end

	def self.memory_length
		@@memory_length
	end

end

class Site_Type

	attr_reader :sites
	attr_reader :title

	@@nall_sites = 0

	def initialize( title, nsites )
		@title = title
		@nsites = nsites
		@bounds = [ @@nall_sites, @@nall_sites+@nsites ]
		@@nall_sites += @nsites
		@sites = []
		init_sites
	end

	def init_site( status )
		@sites << Site.new( status )
	end

	def init_sites
		bounds.each{ |i| init_site( nil ) }
	end

	def update_status( statuses )
		@sites.zip( statuses ){ |site, status| site.update( status ) }
	end

	def number_of_vacancies
		@sites.inject(0){ |count, site| site.is_a_vacancy? ? count += 1 : count }
	end

	def number_occupied
		@sites.inject(0){ |count, site| site.is_a_vacancy? ? count : count += 1 }
	end

	def bounds
		(@bounds[0]...@bounds[1])
	end

	def to_s
		"%8d%8d" % [ number_occupied, number_of_vacancies ]
	end

end

class File

	def get_s
		self.readline.split[0]
	end

	def get_i
		self.readline.split[0].to_i
	end

end

def read_input_from( input_filename )

	# Example input file:
	# sites_atoms.dat
	# 1         nskip
	# 1         number of steps to count for persistent vacancies (set to 0 for no persitence)
	# 3         number of site types
	# tet1      name of site 1
	# 448       number of sites of type 1
	# tet2      name of site type 2
	# 448       number of sites of type 2
	# oct       name of site type 3
	# 448       number of sites of type 3

	abort( "\"#{input_filename}\" not found" ) unless File.exists?( input_filename )

	input_file = File.new( input_filename, 'r' )
	sites_filename = input_file.get_s
	nskip = input_file.get_i
	Site.memory_length = input_file.get_i
	nsite_types = input_file.get_i
	site_types = (0...nsite_types).collect do |i|
		string = input_file.get_s
		nsites = input_file.get_i
		Site_Type.new( string, nsites )
	end

	return sites_filename, nskip, site_types

end

def get_sites_from( string )
	string.split[1..-1].map{ |n| n.to_i }
end

executable_name = File.basename($PROGRAM_NAME)

if ARGV.length == 1
	input_filename = ARGV[0]
else
	abort( "usage: #{executable_name} <input file>" )
end

filename, nskip, site_types = read_input_from( input_filename )
occupancy_data = File.new( filename, 'r' )
step = 0

while( this_line = occupancy_data.gets )
	step += 1
	next unless step%nskip == 0
	data = get_sites_from( this_line )
	site_types.each{ |site_type| site_type.update_status( data[ site_type.bounds ] ) }
	puts "#{step}" + site_types.map{ |site_type| site_type.to_s }.join('   ') if step > Site.memory_length + 1
end
