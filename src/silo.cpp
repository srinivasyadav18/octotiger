/*
 * silo.cpp
 *
 *  Created on: Oct 16, 2018
 *      Author: dmarce1
 */

#define SILO_DRIVER DB_HDF5
#define SILO_VERSION 100

#include "node_registry.hpp"
#include "silo.hpp"
#include "node_server.hpp"
#include "options.hpp"
#include <set>
#include "physcon.hpp"
#include <hpx/lcos/broadcast.hpp>
#include <future>
#include <mutex>
#include "util.hpp"

#define OUTPUT_ROCHE

std::string oct_to_str(int n) {
	char* ptr;
	asprintf( &ptr, "%o", n);
	std::string rc(ptr);
	free(ptr);
	return std::move(rc);
}

struct node_list_t {
	std::vector<node_location::node_id> silo_leaves;
	std::vector<node_location::node_id> all;
	std::vector<integer> positions;
	template<class Arc>
	void serialize(Arc& arc, unsigned) {
		arc & silo_leaves;
		arc & all;
		arc & positions;
	}
};


struct node_entry_t {
	bool load;
	integer position;
	integer locality_id;
	template<class Arc>
	void serialize(Arc& arc, unsigned) {
		arc & load;
		arc & position;
		arc & locality_id;
	}
};


using dir_map_type = std::unordered_map<node_location::node_id, node_entry_t>;
struct mesh_vars_t;
static std::vector<hpx::future<mesh_vars_t>> futs_;
static node_list_t node_list_;
dir_map_type node_dir_;
DBfile* db_;


static const auto& localities = options::all_localities;

static std::mutex silo_mtx_;

void output_stage1(std::string fname, int cycle);
node_list_t output_stage2(std::string fname, int cycle);
void output_stage3(std::string fname, int cycle);

HPX_PLAIN_ACTION(output_stage1, output_stage1_action);
HPX_PLAIN_ACTION(output_stage2, output_stage2_action);
HPX_PLAIN_ACTION(output_stage3, output_stage3_action);

#include <hpx/include/threads.hpp>
#include <hpx/include/run_as.hpp>

template<class T>
struct db_type {
	static constexpr int d = DB_UNKNOWN;
};

template<>
struct db_type<integer> {
	static constexpr int d = DB_LONG_LONG;
};

template<>
struct db_type<real> {
	static constexpr int d = DB_DOUBLE;
};

template<>
struct db_type<char> {
	static constexpr int d = DB_CHAR;
};

template<>
struct db_type<std::int8_t> {
	static constexpr int d = DB_CHAR;
};

constexpr int db_type<integer>::d;
constexpr int db_type<char>::d;
constexpr int db_type<real>::d;

template<class T>
struct write_silo_var {
	void operator()(DBfile* db, const char* name, T var) const {
		int one = 1;
		DBWrite(db, name, &var, &one, 1, db_type<T>::d);
	}
};

template<class T>
struct read_silo_var {
	T operator()(DBfile* db, const char* name) const {
		int one = 1;
		T var;
		if (DBReadVar(db, name, &var) != 0) {
			std::cout << "Unable to read " << name << " \n";
		}
		return var;
	}
};

static const int HOST_NAME_LEN = 100;
static int epoch = 0;
static time_t start_time = time(NULL);
static integer start_step = 0;
static int timestamp;
static int nsteps;
static int time_elapsed;
static int steps_elapsed;
static real output_time;
static real output_rotation_count;

struct mesh_vars_t {
	std::vector<silo_var_t> vars;
	std::vector<std::string> var_names;
	std::vector<std::pair<std::string, real>> outflow;
	std::string mesh_name;
#ifdef OUTPUT_ROCHE
	std::vector<grid::roche_type> roche;
	std::string roche_name;
#endif
	std::vector<std::vector<real>> X;
	std::array<int, NDIM> X_dims;
	std::array<int, NDIM> var_dims;
	node_location location;
	mesh_vars_t(mesh_vars_t&&) = default;
	int compression_level() {
		int lev = 0;
		int i = var_dims[0];
		while (i > INX) {
			i >>= 1;
		}
		return lev;
	}
	mesh_vars_t(const node_location& loc, int compression = 0) :
			location(loc), X(NDIM) {
		const int nx = (INX << compression);
		X_dims[0] = X_dims[1] = X_dims[2] = nx + 1;
		var_dims[0] = var_dims[1] = var_dims[2] = nx;
		const real dx = 2.0 * opts().xscale / nx / (1 << loc.level());
		for (int d = 0; d < NDIM; d++) {
			const int d0 = d;
			X[d0].resize(X_dims[d0]);
			const real o = loc.x_location(d) * opts().xscale;
			for (int i = 0; i < X_dims[d0]; i++) {
				X[d0][i] = (o + i * dx) * opts().code_to_cm;
			}
		}
		mesh_name = oct_to_str(loc.to_id());
#ifdef OUTPUT_ROCHE
		if( opts().problem == DWD) {
			roche.resize(var_dims[0] * var_dims[1] * var_dims[2]);
		}
#endif
	}
};

void output_stage1(std::string fname, int cycle) {
	std::vector<node_location::node_id> ids;
	futs_.clear();
	const auto* node_ptr_ = node_registry::begin()->second.get_ptr().get();
	output_time = node_ptr_->get_time() * opts().code_to_s;
	output_rotation_count = node_ptr_->get_rotation_count();
	for (auto i = node_registry::begin(); i != node_registry::end(); i++) {
		const auto* node_ptr_ = i->second.get_ptr().get();
		if (!node_ptr_->refined()) {
			futs_.push_back(hpx::async([](node_location loc, node_registry::node_ptr ptr)
			{
				const auto* this_ptr = ptr.get_ptr().get();
				assert(this_ptr);
				const real dx = TWO / real(1 << loc.level()) / real(INX);
				mesh_vars_t rc(loc);
				const std::string suffix = oct_to_str(loc.to_id());
				const grid& gridref = this_ptr->get_hydro_grid();
				rc.vars = gridref.var_data();
				rc.outflow = gridref.get_outflows();
#ifdef OUTPUT_ROCHE
				if( opts().problem==DWD ) {
					rc.roche = gridref.get_roche_lobe();
					rc.roche_name = std::string("roche_geometry");
				}
#endif
				return std::move(rc);
			}, i->first, i->second));
		}
	}
}

std::vector<mesh_vars_t> compress(std::vector<mesh_vars_t>&& mesh_vars) {

//	printf("Compressing %i grids...\n", int(mesh_vars.size()));

	using table_type = std::unordered_map<node_location::node_id, std::shared_ptr<mesh_vars_t>>;
	table_type table;
	std::vector<mesh_vars_t> done;
	for (auto& mv : mesh_vars) {
		table.insert(std::make_pair(mv.location.to_id(), std::make_shared<mesh_vars_t>(std::move(mv))));
	}
	int ll = 0;
	std::vector<table_type::iterator> iters;
	while (table.size()) {
		iters.clear();
		auto i = table.begin();
		node_location tmp;
		tmp.from_id(i->first);
		int first_nx = i->second->X_dims[0];
		const integer shift = (tmp.level() - 1) * NDIM;
		if (shift > 0) {
			auto loc_bits = i->first & node_location::node_id(~(0x7 << shift));
			bool levels_match = true;
			for (int j = 0; j < NCHILD; j++) {
				auto this_iter = table.find(loc_bits | (j << shift));
				if (this_iter != table.end()) {
					levels_match = levels_match && (this_iter->second->X_dims[0] == first_nx);
					iters.push_back(this_iter);
				}
			}
			const auto field_names = grid::get_field_names();
			const int nfields = field_names.size();
			if (iters.size() == NCHILD && levels_match) {
				node_location ploc;
				ploc.from_id(loc_bits);
				ploc = ploc.get_parent();
				auto new_mesh_ptr = std::make_shared<mesh_vars_t>(ploc, iters[0]->second->compression_level() + 1);
				for (int f = 0; f < nfields; f++) {
					const bool is_hydro = grid::is_hydro_field(field_names[f]);
					if (is_hydro) {
						new_mesh_ptr->outflow.push_back(std::make_pair(field_names[f], 0.0));
					}
					silo_var_t new_var(field_names[f], new_mesh_ptr->var_dims[0]);
					for (integer ci = 0; ci != NCHILD; ci++) {
						if (is_hydro) {
							new_mesh_ptr->outflow[f].second += iters[ci]->second->outflow[f].second;
						}
						const int nx = new_mesh_ptr->var_dims[0] / 2;
						const int ib = ((ci >> 0) & 1) * nx;
						const int jb = ((ci >> 1) & 1) * nx;
						const int kb = ((ci >> 2) & 1) * nx;
						const int ie = ib + nx;
						const int je = jb + nx;
						const int ke = kb + nx;
						for (int i = ib; i < ie; i++) {
							for (int j = jb; j < je; j++) {
								for (int k = kb; k < ke; k++) {
									const int iiip = i * (4 * nx * nx) + j * 2 * nx + k;
									const int iiic = (i - ib) * (nx * nx) + (j - jb) * nx + (k - kb);
									new_var(iiip) = iters[ci]->second->vars[f](iiic);
								}
							}
						}
					}
					new_mesh_ptr->vars.push_back(std::move(new_var));
				}
#ifdef OUTPUT_ROCHE
				if (opts().problem == DWD) {
					const int nx = new_mesh_ptr->var_dims[0] / 2;
					std::vector<grid::roche_type> roche(8*nx*nx*nx);
					for (integer ci = 0; ci != NCHILD; ci++) {
						const int ib = ((ci >> 0) & 1) * nx;
						const int jb = ((ci >> 1) & 1) * nx;
						const int kb = ((ci >> 2) & 1) * nx;
						const int ie = ib + nx;
						const int je = jb + nx;
						const int ke = kb + nx;
						for (int i = ib; i < ie; i++) {
							for (int j = jb; j < je; j++) {
								for (int k = kb; k < ke; k++) {
									const int iiip = i * (4 * nx * nx) + j * 2 * nx + k;
									const int iiic = (i - ib) * (nx * nx) + (j - jb) * nx + (k - kb);
									roche[iiip] = iters[ci]->second->roche[iiic];
								}
							}
						}
					}
					new_mesh_ptr->roche_name = std::string("roche_geometry");
					new_mesh_ptr->roche = std::move(roche);
				}
#endif
				for (auto this_iter : iters) {
					table.erase(this_iter);
				}
				table.insert(std::make_pair(ploc.to_id(), new_mesh_ptr));
			} else {
				for (auto this_iter : iters) {
					node_location loc;
					loc.from_id(this_iter->first);
					done.push_back(std::move(*this_iter->second));
					table.erase(this_iter);
				}
			}
		} else {
			done.push_back(std::move(*table.begin()->second));
			table.erase(table.begin());
		}
	}

//	printf("Compressed to %i grids...\n", int(done.size()));

	return std::move(done);
}

static std::vector<mesh_vars_t> all_mesh_vars;

node_list_t output_stage2(std::string fname, int cycle) {
	const int this_id = hpx::get_locality_id();
	const int nfields = grid::get_field_names().size();
	std::string this_fname = fname + std::string(".silo");
	all_mesh_vars.clear();
	for (auto& this_fut : futs_) {
		all_mesh_vars.push_back(std::move(GET(this_fut)));
	}
	std::vector<node_location::node_id> ids;
	if (opts().compress_silo) {
		all_mesh_vars = compress(std::move(all_mesh_vars));
	}
	for (const auto& mv : all_mesh_vars) {
		ids.push_back(mv.location.to_id());
	}
	std::vector<node_location::node_id> all;
	std::vector<integer> positions;
	for (auto i = node_registry::begin(); i != node_registry::end(); i++) {
		all.push_back(i->first.to_id());
		positions.push_back(i->second.get_ptr().get()->get_position());
	}
	node_list_t nl;
	nl.silo_leaves = std::move(ids);
	nl.all = std::move(all);
	nl.positions = std::move(positions);
	return std::move(nl);
}

void output_stage3(std::string fname, int cycle) {
	const int this_id = hpx::get_locality_id();
	const int nfields = grid::get_field_names().size();
	std::string this_fname = fname + std::string(".silo");
	double dtime = output_time;
	hpx::threads::run_as_os_thread(
			[&this_fname,this_id,&dtime](integer cycle) {
				DBfile *db;
				if (this_id == 0) {
					db = DBCreateReal(this_fname.c_str(), DB_CLOBBER, DB_LOCAL, "Octo-tiger", SILO_DRIVER);
				} else {
					db = DBOpenReal(this_fname.c_str(), SILO_DRIVER, DB_APPEND);
				}
				float ftime = dtime;
				int one = 1;
				int opt1 = DB_CARTESIAN;
				char cgs[5];
				std::strcpy(cgs,"cgs");
				auto optlist = DBMakeOptlist(9);
				DBAddOption(optlist, DBOPT_HIDE_FROM_GUI, &one);
				DBAddOption(optlist, DBOPT_COORDSYS, &opt1);
				DBAddOption(optlist, DBOPT_CYCLE, &cycle);
				char xstr[2] = {'x', '\0'};
				char ystr[2] = {'y', '\0'};
				char zstr[2] = {'z', '\0'};
				DBAddOption(optlist, DBOPT_XLABEL, xstr );
				DBAddOption(optlist, DBOPT_YLABEL, ystr );
				DBAddOption(optlist, DBOPT_ZLABEL, zstr );
				const char* coord_names[] = {"x", "y", "z"};
				constexpr int data_type = DB_DOUBLE;
				constexpr int ndim = NDIM;
				constexpr int coord_type = DB_COLLINEAR;
				int count = 0;
				for( const auto& mesh_vars : all_mesh_vars) {
					const auto& X = mesh_vars.X;
					const real* coords[NDIM];
					for( int d = 0; d < NDIM; d++) {
						coords[d] = X[d].data();
					}
					const auto& dir_name = mesh_vars.mesh_name.c_str();
					DBMkDir(db, dir_name);
					DBSetDir( db, dir_name );
					DBPutQuadmesh(db, "quadmesh", coord_names, coords, mesh_vars.X_dims.data(), ndim, data_type, coord_type, optlist);
					for ( integer m = 0; m != mesh_vars.vars.size(); m++) {
						const auto& o = mesh_vars.vars[m];
						const bool is_hydro = grid::is_hydro_field(o.name());
						if( is_hydro ) {
							real outflow = mesh_vars.outflow[m].second;
							DBAddOption(optlist, DBOPT_DTIME, &outflow);
							DBAddOption(optlist, DBOPT_CONSERVED, &one);
						}
						DBPutQuadvar1(db, o.name(), "quadmesh", o.data(), mesh_vars.var_dims.data(), ndim, (const void*) NULL, 0,
								DB_DOUBLE, DB_ZONECENT, optlist);
						count++;
						if( is_hydro ) {
							DBClearOption(optlist, DBOPT_CONSERVED);
							DBClearOption(optlist, DBOPT_DTIME);
						}
					}
#ifdef OUTPUT_ROCHE
					if( opts().problem==DWD) {
						auto this_name = mesh_vars.roche_name;
						DBPutQuadvar1(db, this_name.c_str(), mesh_vars.mesh_name.c_str(), mesh_vars.roche.data(), mesh_vars.var_dims.data(), ndim, (const void*) NULL, 0,
								db_type<grid::roche_type>::d, DB_ZONECENT, optlist);
					}
#endif
					DBSetDir( db, "/" );
				}
				DBFreeOptlist( optlist);
				DBClose( db);
			}, cycle).get();
	if (this_id < integer(localities.size()) - 1) {
		output_stage3_action func;
		func(localities[this_id + 1], fname, cycle);
	}

	double rtime = output_rotation_count;
	if (this_id == 0) {
		hpx::threads::run_as_os_thread(
				[&this_fname,nfields,&rtime](int cycle) {
					auto* db = DBOpenReal(this_fname.c_str(), SILO_DRIVER, DB_APPEND);
					double dtime = output_time;
					float ftime = dtime;
					std::vector<node_location> node_locs;
					std::vector<char*> mesh_names;
					std::vector<std::vector<char*>> field_names(nfields);
#ifdef OUTPUT_ROCHE
					std::vector<char*> roche_names;
#endif
					for (auto& i : node_list_.silo_leaves) {
						node_location nloc;
						nloc.from_id(i);
						node_locs.push_back(nloc);
					}
					const auto top_field_names = grid::get_field_names();
					for (int i = 0; i < node_locs.size(); i++) {
						const auto suffix = oct_to_str(node_locs[i].to_id());
						const auto str = "/" + suffix + "/quadmesh";
						char* ptr = new char[str.size() + 1];
						std::strcpy(ptr, str.c_str());
						mesh_names.push_back(ptr);
						for (int f = 0; f < nfields; f++) {
							const auto str = "/" + suffix + "/" + top_field_names[f];
							char* ptr = new char[str.size() + 1];
							strcpy(ptr, str.c_str());
							field_names[f].push_back(ptr);
						}
#ifdef OUTPUT_ROCHE
						if( opts().problem == DWD ) {
							const auto str = "/" + suffix + std::string("/roche_geometry");
							char* ptr = new char[str.size() + 1];
							strcpy(ptr, str.c_str());
							roche_names.push_back(ptr);
						}
#endif
					}

					const int n_total_domains = mesh_names.size();

					auto optlist = DBMakeOptlist(9);
					int opt1 = DB_CARTESIAN;
					int mesh_type = DB_QUADMESH;
					char cgs[5];
					std::strcpy(cgs,"cgs");
					DBAddOption(optlist, DBOPT_COORDSYS, &opt1);
					DBAddOption(optlist, DBOPT_CYCLE, &cycle);
					DBAddOption(optlist, DBOPT_DTIME, &dtime);
					DBAddOption(optlist, DBOPT_TIME, &ftime);
					DBAddOption(optlist, DBOPT_MB_BLOCK_TYPE, &mesh_type);
					DBAddOption(optlist, DBOPT_XUNITS, cgs);
					char xstr[2] = {'x', '\0'};
					char ystr[2] = {'y', '\0'};
					char zstr[2] = {'z', '\0'};
					DBAddOption(optlist, DBOPT_XLABEL, xstr );
					DBAddOption(optlist, DBOPT_YLABEL, ystr );
					DBAddOption(optlist, DBOPT_ZLABEL, zstr );
					assert( n_total_domains > 0 );
					printf( "Putting %i\n", n_total_domains );
					DBPutMultimesh(db, "quadmesh", n_total_domains, mesh_names.data(), NULL, optlist);
					for (int f = 0; f < nfields; f++) {
						DBPutMultivar( db, top_field_names[f].c_str(), n_total_domains, field_names[f].data(), std::vector<int>(n_total_domains, DB_QUADVAR).data(), optlist);
					}
#ifdef OUTPUT_ROCHE
					if( opts().problem == DWD ) {
						DBPutMultivar( db, "roche_geometry", n_total_domains, roche_names.data(), std::vector<int>(n_total_domains, DB_QUADVAR).data(), optlist);
					}
#endif
					write_silo_var<integer> fi;
					write_silo_var<real> fr;
					fi(db, "version", SILO_VERSION);
					fr(db, "code_to_g", opts().code_to_g);
					fr(db, "code_to_s", opts().code_to_s);
					fr(db, "code_to_cm", opts().code_to_cm);
					fi(db, "n_species", integer(opts().n_species));
					fi(db, "eos", integer(opts().eos));
					fi(db, "gravity", integer(opts().gravity));
					fi(db, "hydro", integer(opts().hydro));
					fr(db, "omega", grid::get_omega() / opts().code_to_s);
					fr(db, "output_frequency", opts().output_dt);
					fi(db, "problem", integer(opts().problem));
					fi(db, "radiation", integer(opts().radiation));
					fr(db, "refinement_floor", opts().refinement_floor);
					fr(db, "cgs_time", dtime);
					fr(db, "rotational_time", rtime);
					fr(db, "xscale", opts().xscale); char hostname[HOST_NAME_LEN];
					gethostname(hostname, HOST_NAME_LEN);
					DBWrite( db, "hostname", hostname, &HOST_NAME_LEN, 1, DB_CHAR);
					int nnodes = node_list_.all.size();
					DBWrite( db, "node_list", node_list_.all.data(), &nnodes, 1, DB_LONG_LONG);
					DBWrite( db, "node_positions", node_list_.positions.data(), &nnodes, 1, db_type<integer>::d);
					int nspc = opts().n_species;
					DBWrite( db, "X", opts().X.data(), &nspc, 1, db_type<real>::d);
					DBWrite( db, "Z", opts().Z.data(), &nspc, 1, db_type<real>::d);
					DBWrite( db, "atomic_mass", opts().atomic_mass.data(), &nspc, 1, db_type<real>::d);
					DBWrite( db, "atomic_number", opts().atomic_number.data(), &nspc, 1, db_type<real>::d);
					fi(db, "node_count", integer(nnodes));
					fi(db, "leaf_count", integer(node_list_.silo_leaves.size()));
					write_silo_var<integer>()(db, "timestamp", timestamp);
					write_silo_var<integer>()(db, "epoch", epoch);
					write_silo_var<integer>()(db, "locality_count", localities.size());
					write_silo_var<integer>()(db, "thread_count", localities.size() * std::thread::hardware_concurrency());
					write_silo_var<integer>()(db, "step_count", nsteps);
					write_silo_var<integer>()(db, "time_elapsed", time_elapsed);
					write_silo_var<integer>()(db, "steps_elapsed", steps_elapsed);

					DBFreeOptlist( optlist);
					DBClose( db);
					for (auto ptr : mesh_names) {
						delete[] ptr;
					}
#ifdef OUTPUT_ROCHE
					if( opts().problem == DWD ) {
						for (auto ptr : roche_names) {
							delete[] ptr;
						}
					}
#endif
					for (int f = 0; f < nfields; f++) {
						for (auto& s : field_names[f]) {
							delete[] s;
						}
					}
				}, cycle).get();
	}
}

void output_all(std::string fname, int cycle, bool block) {

	static hpx::future<void> barrier(hpx::make_ready_future<void>());
	GET(barrier);
	nsteps = node_registry::begin()->second.get_ptr().get()->get_step_num();
	timestamp = time(NULL);
	steps_elapsed = nsteps - start_step;
	time_elapsed = time(NULL) - start_time;
	start_time = timestamp;
	start_step = nsteps;
	std::vector<hpx::future<void>> futs1;
	for (auto& id : localities) {
		futs1.push_back(hpx::async<output_stage1_action>(id, fname, cycle));
	}
	GET(hpx::when_all(futs1));

	std::vector<hpx::future<node_list_t>> id_futs;
	for (auto& id : localities) {
		id_futs.push_back(hpx::async<output_stage2_action>(id, fname, cycle));
	}
	node_list_.silo_leaves.clear();
	node_list_.all.clear();
	node_list_.positions.clear();
	for (auto& f : id_futs) {
		node_list_t this_list = GET(f);
		for (auto& i : this_list.silo_leaves) {
			node_list_.silo_leaves.push_back(i);
		}
		for (auto& i : this_list.all) {
			node_list_.all.push_back(i);
		}
		for (auto& i : this_list.positions) {
			node_list_.positions.push_back(i);
		}
	}
    barrier = hpx::async([cycle,fname]() {
		GET(hpx::async<output_stage3_action>(localities[0], fname, cycle));
	});

	if (block) {
		GET(barrier);
		barrier = hpx::make_ready_future<void>();
	}
}

#define SILO_TEST(i) \
	if( i != 0 ) printf( "SILO call failed at %i\n", __LINE__ );

void load_options_from_silo(std::string fname, DBfile* db) {
	const auto func =
			[&fname,&db]()
			{
				bool leaveopen;
				if( db == NULL )
				{
					db = DBOpenReal( fname.c_str(), SILO_DRIVER, DB_READ);
					leaveopen = false;
				}
				else
				{
					leaveopen = true;
				}
				if (db != NULL)
				{

					read_silo_var<integer> ri;
					read_silo_var<real> rr;
					integer version = ri( db, "version");
					if( version > SILO_VERSION) {
						printf( "WARNING: Attempting to load a version %i SILO file, maximum version allowed for this Octo-tiger is %i\n", int(version), SILO_VERSION);
					}
					opts().code_to_g = rr(db, "code_to_g");
					opts().code_to_s = rr(db, "code_to_s");
					opts().code_to_cm = rr(db, "code_to_cm");
					opts().n_species = ri(db, "n_species");
					opts().eos = eos_type(ri(db, "eos"));
					opts().gravity = ri(db, "gravity");
					opts().hydro = ri(db, "hydro");
					opts().omega = rr(db, "omega") * opts().code_to_s;
					opts().output_dt = rr(db, "output_frequency");
					opts().problem = problem_type(ri(db, "problem"));
					opts().radiation = ri(db, "radiation");
					opts().refinement_floor = rr(db, "refinement_floor");
					opts().xscale = rr(db, "xscale");
					opts().atomic_number.resize(opts().n_species);
					opts().atomic_mass.resize(opts().n_species);
					opts().X.resize(opts().n_species);
					opts().Z.resize(opts().n_species);
					SILO_TEST(DBReadVar(db, "atomic_number", opts().atomic_number.data()));
					SILO_TEST(DBReadVar(db, "atomic_mass", opts().atomic_mass.data()));
					SILO_TEST(DBReadVar(db, "X", opts().X.data()));
					SILO_TEST(DBReadVar(db, "Z", opts().Z.data()));
					if( !leaveopen)
					{
						DBClose(db);
					}
				}
				else
				{
					std::cout << "Could not load " << fname;
					throw;
				}
			};
	if (db == NULL) {
		GET(hpx::threads::run_as_os_thread(func));
	} else {
		func();
	}
	grid::set_omega(opts().omega, false);
	set_units(opts().code_to_g, opts().code_to_cm, opts().code_to_s, 1); /**/

}

void load_open(std::string fname, dir_map_type map) {
	printf( "LOAD OPENED on proc %i\n", hpx::get_locality_id());
	load_options_from_silo(fname, db_); /**/
	hpx::threads::run_as_os_thread([&]() {
		db_ = DBOpenReal( fname.c_str(), SILO_DRIVER, DB_READ);
		read_silo_var<real> rr;
		output_time = rr(db_, "cgs_time"); /**/
		output_rotation_count = 2 * M_PI * rr(db_, "rotational_time"); /**/
		printf( "rotational_time = %e\n", output_rotation_count);
		output_time /= opts().code_to_s;
		node_dir_ = std::move(map);
		printf( "%e\n", output_time );
//		sleep(100);
	}).get();
}

void load_close() {
	DBClose(db_);
}

HPX_PLAIN_ACTION(load_close, load_close_action);
HPX_PLAIN_ACTION(load_open, load_open_action);

node_server::node_server(const node_location& loc) :
		my_location(loc) {
	const auto& localities = opts().all_localities;
	initialize(0.0, 0.0);
	step_num = gcycle = hcycle = rcycle = 0;

	auto iter = node_dir_.find(loc.to_id());
	assert(iter != node_dir_.end());

	if (!iter->second.load) {
		printf("Creating %s on %i\n", loc.to_str().c_str(), int(hpx::get_locality_id()));
		int nc = 0;
		for (int ci = 0; ci < NCHILD; ci++) {
			auto cloc = loc.get_child(ci);
			auto iter = node_dir_.find(cloc.to_id());
			if (iter != node_dir_.end() && iter->second.locality_id != hpx::get_locality_id()) {
				is_refined = true;
				children[ci] = hpx::new_<node_server>(localities[iter->second.locality_id], cloc);
				nc++;
			}
		}
		for (int ci = 0; ci < NCHILD; ci++) {
			auto cloc = loc.get_child(ci);
			auto iter = node_dir_.find(cloc.to_id());
			if (iter != node_dir_.end() && iter->second.locality_id == hpx::get_locality_id()) {
				is_refined = true;
				children[ci] = hpx::new_<node_server>(localities[iter->second.locality_id], cloc);
				nc++;
			}
		}
		assert(nc == 0 || nc == NCHILD);
	} else {
		printf("Loading %s on %i\n", loc.to_str().c_str(), int(hpx::get_locality_id()));
		silo_load_t load;
		static const auto hydro_names = grid::get_hydro_field_names();
		load.vars.resize(hydro_names.size());
		load.outflows.resize(hydro_names.size());
		hpx::threads::run_as_os_thread([&]() {
			static std::mutex mtx;
			std::lock_guard<std::mutex> lock(mtx);
			const std::string suffix = oct_to_str(loc.to_id());
			for( int f = 0; f != hydro_names.size(); f++) {
				const auto this_name = std::string( "/") + suffix + std::string( "/") + hydro_names[f]; /**/
				auto var = DBGetQuadvar(db_,this_name.c_str());
				load.nx = var->dims[0];
				const int nvar = load.nx * load.nx * load.nx;
				load.outflows[f].first = load.vars[f].first = hydro_names[f];
				load.vars[f].second.resize(nvar);
				load.outflows[f].second = var->dtime;
				std::memcpy(load.vars[f].second.data(), var->vals[0], sizeof(real)*nvar);
				DBFreeQuadvar(var);
			}
		}).get();
		if (load.nx == INX) {
			is_refined = false;
			for (integer f = 0; f < hydro_names.size(); f++) {
				grid_ptr->set(load.vars[f].first, load.vars[f].second.data());
				grid_ptr->set_outflow(std::move(load.outflows[f]));
			}
			grid_ptr->rho_from_species();
		} else {
			is_refined = true;
			auto child_loads = load.decompress();
			for (integer ci = 0; ci < NCHILD; ci++) {
				auto cloc = loc.get_child(ci);
				auto iter = node_dir_.find(cloc.to_id());
				children[ci] = hpx::new_<node_server>(localities[iter->second.locality_id], cloc, child_loads[ci]);
			}
		}
	}
	current_time = output_time;
	rotational_time = output_rotation_count;
}

node_server::node_server(const node_location& loc, silo_load_t load) :
		my_location(loc) {
	printf("Distributing %s on %i\n", loc.to_str().c_str(), int(hpx::get_locality_id()));
	const auto& localities = opts().all_localities;
	initialize(0.0, 0.0);
	step_num = gcycle = hcycle = rcycle = 0;
	int nc = 0;
	static const auto hydro_names = grid::get_hydro_field_names();
	if (load.nx == INX) {
		is_refined = false;
		for (integer f = 0; f < hydro_names.size(); f++) {
			grid_ptr->set(load.vars[f].first, load.vars[f].second.data());
			grid_ptr->set_outflow(std::move(load.outflows[f]));
		}
		grid_ptr->rho_from_species();
	} else {
		is_refined = true;
		auto child_loads = load.decompress();
		for (integer ci = 0; ci < NCHILD; ci++) {
			auto cloc = loc.get_child(ci);
			auto iter = node_dir_.find(cloc.to_id());
			children[ci] = hpx::new_<node_server>(localities[iter->second.locality_id], cloc, child_loads[ci]);
		}
	}
	current_time = output_time;
	rotational_time = output_rotation_count;
	assert(nc == 0 || nc == NCHILD);
}

void load_data_from_silo(std::string fname, node_server* root_ptr, hpx::id_type root) {
	const integer nprocs = opts().all_localities.size();
	static int sz = localities.size();
	DBfile* db = GET(hpx::threads::run_as_os_thread(DBOpenReal, fname.c_str(), SILO_DRIVER, DB_READ));
	epoch = GET(hpx::threads::run_as_os_thread(read_silo_var<integer>(), db, "epoch"));
	epoch++;
	std::vector<node_location::node_id> node_list;
	std::vector<integer> positions;
	std::vector<hpx::future<void>> futs;
	int node_count;
	if (db != NULL) {
		DBmultimesh* master_mesh = GET(hpx::threads::run_as_os_thread([&]()
		{
			return DBGetMultimesh( db, "quadmesh");
		}));
		const int chunk_size = std::ceil(real(master_mesh->nblocks) / real(sz));
		hpx::threads::run_as_os_thread([&]() {
			const read_silo_var<integer> ri;
			node_count = ri(db,"node_count");
			node_list.resize(node_count);
			positions.resize(node_count);
			DBReadVar(db, "node_list", node_list.data());
			DBReadVar(db, "node_positions", positions.data());
		}).get();
		GET(hpx::threads::run_as_os_thread(DBClose, db));
		std::set<node_location::node_id> load_locs;
		for (int i = 0; i < master_mesh->nblocks; i++) {
			const int num = std::strtol(master_mesh->meshnames[i] + 1, nullptr, 8);
			printf("%i\n", num);
			load_locs.insert(num);
		}
		for (int i = 0; i < node_list.size(); i++) {
			node_entry_t entry;
			entry.position = positions[i];
			entry.load = bool(load_locs.find(node_list[i]) != load_locs.end());
			entry.locality_id = positions[i] * nprocs / positions.size();
			node_dir_[node_list[i]] = entry;
		}
		auto this_dir = std::move(node_dir_);
		for (int i = 0; i < nprocs; i++) {
			printf( "Sending LOAD OPEN to %i\n", i);
			futs.push_back(hpx::async<load_open_action>(opts().all_localities[i], fname, this_dir));
		}
		GET(hpx::threads::run_as_os_thread(DBFreeMultimesh, master_mesh));
		for (auto& f : futs) {
			GET(f);
		}
	} else {
		std::cout << "Could not load " << fname;
		throw;
	}
	root_ptr->reconstruct_tree();
	node_registry::clear();
	futs.clear();
	for (int i = 0; i < nprocs; i++) {
		futs.push_back(hpx::async<load_close_action>(opts().all_localities[i]));
	}
	for (auto& f : futs) {
		GET(f);
	}
}

void node_server::reconstruct_tree() {
	for (integer ci = 0; ci < NCHILD; ci++) {
		is_refined = true;
		auto cloc = my_location.get_child(ci);
		auto iter = node_dir_.find(cloc.to_id());
		children[ci] = hpx::new_<node_server>(localities[iter->second.locality_id], cloc).get();
	}
	current_time = output_time;
	rotational_time = output_rotation_count;
}

silo_var_t::silo_var_t(const std::string& name, std::size_t nx) :
		name_(name), data_(nx * nx * nx) {
}

double&
silo_var_t::operator()(int i) {
	return data_[i];
}

double silo_var_t::operator()(int i) const {
	return data_[i];
}

std::vector<silo_load_t> silo_load_t::decompress() {
	std::vector<silo_load_t> children;
	assert(nx > INX);
	for (int ci = 0; ci < NCHILD; ci++) {
		silo_load_t child;
		child.nx = nx / 2;
		child.vars.resize(vars.size());
		child.outflows.resize(vars.size());
		const integer xo = (ci & (1 << XDIM)) ? child.nx : 0;
		const integer yo = (ci & (1 << YDIM)) ? child.nx : 0;
		const integer zo = (ci & (1 << ZDIM)) ? child.nx : 0;
		for (int f = 0; f < vars.size(); f++) {
			child.vars[f].second.resize(child.nx * child.nx * child.nx);
			child.outflows[f].first = child.vars[f].first = vars[f].first;
			child.outflows[f].second = outflows[f].second / NCHILD;
			for (integer cx = 0; cx < child.nx; cx++) {
				for (integer cy = 0; cy < child.nx; cy++) {
					for (integer cz = 0; cz < child.nx; cz++) {
						const integer child_index = cx + child.nx * (cy + child.nx * cz);
						const integer parent_index = (cx + xo) + nx * ((cy + yo) + nx * (cz + zo));
						child.vars[f].second[child_index] = vars[f].second[parent_index];
					}
				}
			}
		}

		children.push_back(std::move(child));
	}
	vars.clear();
	outflows.clear();
	return std::move(children);
}


