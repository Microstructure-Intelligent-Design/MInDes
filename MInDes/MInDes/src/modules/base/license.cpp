#include "license.h"

#include "MACRO_DEF.h"
#include "license_public_key.h"
#include "whereami.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <system_error>

#if defined(_WIN32)
#define NOMINMAX
#include <windows.h>
#include <bcrypt.h>
#include <ncrypt.h>
#pragma comment(lib, "bcrypt.lib")
#pragma comment(lib, "ncrypt.lib")
#elif defined(__linux__)
#include <openssl/bio.h>
#include <openssl/ec.h>
#include <openssl/ecdsa.h>
#include <openssl/evp.h>
#include <openssl/pem.h>
#include <openssl/rand.h>
#include <openssl/rsa.h>
#if defined(MINDES_HAS_TPM2_TSS)
#include <tss2/tss2_fapi.h>
#endif
#else
#error "MInDes License V2 supports only Windows and Linux"
#endif

namespace {
	using Bytes = std::vector<std::uint8_t>;
	using Hash256 = std::array<std::uint8_t, 32>;
	using Id128 = std::array<std::uint8_t, 16>;

	constexpr std::uint16_t format_version = 2;
	constexpr std::uint16_t binding_version = 2;
	constexpr std::size_t magic_size = 20;
	constexpr std::size_t maximum_file_size = 1024 * 1024;
	constexpr char legacy_magic[] = "## MINDES PROGRAM LICENSE ##";
	constexpr char active_magic[] = "MINDES-ACTIVE-V2";
	constexpr char request_magic[] = "MINDES-REQUEST-V2";
	constexpr char license_magic[] = "MINDES-LICENSE-V2";
	constexpr char tpm_domain[] = "MINDES-TPM-BINDING-V2";
	constexpr char smbios_domain[] = "MINDES-SMBIOS-BINDING-V2";

	class Writer {
	public:
		void bytes(const void* data, std::size_t size) {
			const auto* first = static_cast<const std::uint8_t*>(data);
			data_.insert(data_.end(), first, first + size);
		}

		template <std::size_t N>
		void bytes(const std::array<std::uint8_t, N>& value) {
			bytes(value.data(), value.size());
		}

		void vector(const Bytes& value) {
			u32(static_cast<std::uint32_t>(value.size()));
			bytes(value.data(), value.size());
		}

		void string(const std::string& value) {
			u32(static_cast<std::uint32_t>(value.size()));
			bytes(value.data(), value.size());
		}

		void u8(std::uint8_t value) { data_.push_back(value); }
		void u16(std::uint16_t value) {
			for (int shift = 0; shift < 16; shift += 8) u8(static_cast<std::uint8_t>(value >> shift));
		}
		void u32(std::uint32_t value) {
			for (int shift = 0; shift < 32; shift += 8) u8(static_cast<std::uint8_t>(value >> shift));
		}
		void u64(std::uint64_t value) {
			for (int shift = 0; shift < 64; shift += 8) u8(static_cast<std::uint8_t>(value >> shift));
		}
		void i64(std::int64_t value) { u64(static_cast<std::uint64_t>(value)); }

		void magic(const char* value) {
			std::array<std::uint8_t, magic_size> padded{};
			const std::size_t length = std::min(std::strlen(value), padded.size());
			std::memcpy(padded.data(), value, length);
			bytes(padded);
		}

		const Bytes& data() const { return data_; }
		Bytes take() { return std::move(data_); }

	private:
		Bytes data_;
	};

	class Reader {
	public:
		explicit Reader(const Bytes& data) : data_(data) {}

		bool bytes(void* destination, std::size_t size) {
			if (size > data_.size() - offset_) return false;
			std::memcpy(destination, data_.data() + offset_, size);
			offset_ += size;
			return true;
		}

		template <std::size_t N>
		bool bytes(std::array<std::uint8_t, N>& value) { return bytes(value.data(), value.size()); }

		bool vector(Bytes& value, std::size_t maximum = 65536) {
			std::uint32_t size = 0;
			if (!u32(size) || size > maximum || size > data_.size() - offset_) return false;
			value.assign(data_.begin() + static_cast<std::ptrdiff_t>(offset_),
				data_.begin() + static_cast<std::ptrdiff_t>(offset_ + size));
			offset_ += size;
			return true;
		}

		bool string(std::string& value, std::size_t maximum = 4096) {
			Bytes data;
			if (!vector(data, maximum)) return false;
			value.assign(data.begin(), data.end());
			return true;
		}

		bool u8(std::uint8_t& value) {
			if (offset_ >= data_.size()) return false;
			value = data_[offset_++];
			return true;
		}
		bool u16(std::uint16_t& value) {
			value = 0;
			for (int shift = 0; shift < 16; shift += 8) { std::uint8_t byte = 0; if (!u8(byte)) return false; value |= static_cast<std::uint16_t>(byte) << shift; }
			return true;
		}
		bool u32(std::uint32_t& value) {
			value = 0;
			for (int shift = 0; shift < 32; shift += 8) { std::uint8_t byte = 0; if (!u8(byte)) return false; value |= static_cast<std::uint32_t>(byte) << shift; }
			return true;
		}
		bool u64(std::uint64_t& value) {
			value = 0;
			for (int shift = 0; shift < 64; shift += 8) { std::uint8_t byte = 0; if (!u8(byte)) return false; value |= static_cast<std::uint64_t>(byte) << shift; }
			return true;
		}
		bool i64(std::int64_t& value) { std::uint64_t raw = 0; if (!u64(raw)) return false; value = static_cast<std::int64_t>(raw); return true; }

		bool magic(const char* expected) {
			std::array<std::uint8_t, magic_size> actual{};
			if (!bytes(actual)) return false;
			std::array<std::uint8_t, magic_size> wanted{};
			std::memcpy(wanted.data(), expected, std::min(std::strlen(expected), wanted.size()));
			return actual == wanted;
		}

		std::size_t offset() const { return offset_; }
		bool finished() const { return offset_ == data_.size(); }

	private:
		const Bytes& data_;
		std::size_t offset_ = 0;
	};

	struct CommonRecord {
		pf::BindingType binding_type = pf::BindingType::Unknown;
		pf::DeviceKeyAlgorithm key_algorithm = pf::DeviceKeyAlgorithm::None;
		Hash256 machine_hash{};
		Id128 request_id{};
	};

	struct ActiveRecord {
		CommonRecord common;
		Hash256 request_digest{};
		std::int64_t created_at = 0;
		std::int64_t last_seen = 0;
		std::uint64_t sequence = 0;
		Bytes public_key;
		std::string key_reference;
		std::string fallback_reason;
		Bytes signature;
		Bytes signed_body;
	};

	struct RequestRecord {
		CommonRecord common;
		Hash256 nonce{};
		std::int64_t created_at = 0;
		Bytes public_key;
		std::string fallback_reason;
		Bytes proof;
		Bytes signed_body;
	};

	struct LicenseRecord {
		CommonRecord common;
		Hash256 request_digest{};
		Id128 license_id{};
		std::int64_t issued_at = 0;
		std::int64_t expires_at = 0;
		std::uint64_t features = 0;
		Bytes public_key;
		std::uint8_t issuer_algorithm = 1;
		Bytes signature;
		Bytes signed_body;
	};

	struct BindingMaterial {
		pf::MachineIdentity identity;
		std::string key_reference;
		std::string fallback_reason;
	};

	std::int64_t unix_time_now() {
		return std::chrono::duration_cast<std::chrono::seconds>(
			std::chrono::system_clock::now().time_since_epoch()).count();
	}

	std::filesystem::path get_executable_directory() {
		const int required = wai_getExecutablePath(nullptr, 0, nullptr);
		if (required <= 0) throw std::runtime_error("Failed to get executable path length");
		std::vector<char> path(static_cast<std::size_t>(required) + 1, '\0');
		const int written = wai_getExecutablePath(path.data(), required, nullptr);
		if (written <= 0 || written > required) throw std::runtime_error("Failed to get executable path");
		path[static_cast<std::size_t>(written)] = '\0';
		std::error_code error;
		auto executable = std::filesystem::canonical(std::filesystem::u8path(path.data()), error);
		if (error) throw std::runtime_error("Failed to canonicalize executable path: " + error.message());
		auto directory = executable.parent_path();
		if (directory.empty() || !std::filesystem::is_directory(directory, error) || error)
			throw std::runtime_error("Invalid executable directory: " + directory.string());
		return directory;
	}

	bool read_file(const std::filesystem::path& path, Bytes& data) {
		std::ifstream input(path, std::ios::binary | std::ios::ate);
		if (!input) return false;
		const auto size = input.tellg();
		if (size < 0 || static_cast<std::uint64_t>(size) > maximum_file_size) return false;
		data.resize(static_cast<std::size_t>(size));
		input.seekg(0);
		return data.empty() || static_cast<bool>(input.read(
			reinterpret_cast<char*>(data.data()), static_cast<std::streamsize>(size)));
	}

	bool random_bytes(std::uint8_t* output, std::size_t size);

	bool atomic_write(const std::filesystem::path& path, const Bytes& data, std::string& error_message) {
		std::array<std::uint8_t, 8> temporary_id{};
		if (!random_bytes(temporary_id.data(), temporary_id.size())) {
			error_message = "Unable to create a temporary filename for: " + path.string();
			return false;
		}
		static constexpr char hex[] = "0123456789ABCDEF";
		std::string suffix = ".tmp.";
		for (const auto value : temporary_id) {
			suffix.push_back(hex[value >> 4]);
			suffix.push_back(hex[value & 0x0F]);
		}
		std::filesystem::path temporary = path;
		temporary += suffix;
		{
			std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
			if (!output || (!data.empty() && !output.write(reinterpret_cast<const char*>(data.data()), static_cast<std::streamsize>(data.size())))) {
				error_message = "Unable to write temporary file: " + temporary.string();
				return false;
			}
			output.flush();
			if (!output) { error_message = "Unable to flush temporary file: " + temporary.string(); return false; }
		}
#if defined(_WIN32)
		if (!MoveFileExW(temporary.c_str(), path.c_str(), MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH)) {
			error_message = "Unable to replace file atomically: " + path.string();
			DeleteFileW(temporary.c_str());
			return false;
		}
#else
		std::error_code error;
		std::filesystem::rename(temporary, path, error);
		if (error) { error_message = "Unable to replace file atomically: " + path.string() + ": " + error.message(); std::filesystem::remove(temporary); return false; }
#endif
		return true;
	}

	bool random_bytes(std::uint8_t* output, std::size_t size) {
#if defined(_WIN32)
		return BCryptGenRandom(nullptr, output, static_cast<ULONG>(size), BCRYPT_USE_SYSTEM_PREFERRED_RNG) == 0;
#else
		return RAND_bytes(output, static_cast<int>(size)) == 1;
#endif
	}

	bool sha256(const std::uint8_t* data, std::size_t size, Hash256& digest) {
#if defined(_WIN32)
		BCRYPT_ALG_HANDLE algorithm = nullptr;
		BCRYPT_HASH_HANDLE hash = nullptr;
		DWORD object_size = 0, returned = 0;
		Bytes object;
		bool success = BCryptOpenAlgorithmProvider(&algorithm, BCRYPT_SHA256_ALGORITHM, nullptr, 0) == 0
			&& BCryptGetProperty(algorithm, BCRYPT_OBJECT_LENGTH, reinterpret_cast<PUCHAR>(&object_size), sizeof(object_size), &returned, 0) == 0;
		if (success) object.resize(object_size);
		success = success && BCryptCreateHash(algorithm, &hash, object.data(), object_size, nullptr, 0, 0) == 0
			&& BCryptHashData(hash, const_cast<PUCHAR>(data), static_cast<ULONG>(size), 0) == 0
			&& BCryptFinishHash(hash, digest.data(), static_cast<ULONG>(digest.size()), 0) == 0;
		if (hash) BCryptDestroyHash(hash);
		if (algorithm) BCryptCloseAlgorithmProvider(algorithm, 0);
		return success;
#else
		unsigned int length = 0;
		return EVP_Digest(data, size, digest.data(), &length, EVP_sha256(), nullptr) == 1 && length == digest.size();
#endif
	}

	bool sha256(const Bytes& data, Hash256& digest) { return sha256(data.data(), data.size(), digest); }

	std::string hex_string(const std::uint8_t* data, std::size_t size) {
		std::ostringstream output;
		output << std::uppercase << std::hex << std::setfill('0');
		for (std::size_t index = 0; index < size; ++index) output << std::setw(2) << static_cast<unsigned int>(data[index]);
		return output.str();
	}

	std::string display_code(const Hash256& hash) { return hex_string(hash.data(), 16); }

	void write_common(Writer& writer, const char* magic, const CommonRecord& record) {
		writer.magic(magic);
		writer.u16(format_version);
		writer.u16(binding_version);
		writer.u8(static_cast<std::uint8_t>(record.binding_type));
		writer.u8(static_cast<std::uint8_t>(record.key_algorithm));
		writer.u16(0);
		writer.bytes(record.machine_hash);
		writer.bytes(record.request_id);
	}

	bool read_common(Reader& reader, const char* magic, CommonRecord& record) {
		std::uint16_t file_version = 0, bind_version = 0, reserved = 0;
		std::uint8_t binding = 0, algorithm = 0;
		if (!reader.magic(magic) || !reader.u16(file_version) || !reader.u16(bind_version)
			|| !reader.u8(binding) || !reader.u8(algorithm) || !reader.u16(reserved)
			|| !reader.bytes(record.machine_hash) || !reader.bytes(record.request_id)) return false;
		if (file_version != format_version || bind_version != binding_version || reserved != 0) return false;
		if (binding != static_cast<std::uint8_t>(pf::BindingType::Tpm2Key)
			&& binding != static_cast<std::uint8_t>(pf::BindingType::Smbios)) return false;
		record.binding_type = static_cast<pf::BindingType>(binding);
		record.key_algorithm = static_cast<pf::DeviceKeyAlgorithm>(algorithm);
		return true;
	}

	Bytes serialize_request_body(const RequestRecord& value) {
		Writer writer;
		write_common(writer, request_magic, value.common);
		writer.bytes(value.nonce);
		writer.i64(value.created_at);
		writer.vector(value.public_key);
		writer.string(value.fallback_reason);
		return writer.take();
	}

	Bytes serialize_request(const RequestRecord& value) {
		Writer writer;
		const Bytes body = serialize_request_body(value);
		writer.bytes(body.data(), body.size());
		writer.vector(value.proof);
		return writer.take();
	}

	bool parse_request(const Bytes& data, RequestRecord& value) {
		Reader reader(data);
		if (!read_common(reader, request_magic, value.common) || !reader.bytes(value.nonce)
			|| !reader.i64(value.created_at) || !reader.vector(value.public_key)
			|| !reader.string(value.fallback_reason)) return false;
		value.signed_body.assign(data.begin(), data.begin() + static_cast<std::ptrdiff_t>(reader.offset()));
		return reader.vector(value.proof) && reader.finished();
	}

	Bytes serialize_active_body(const ActiveRecord& value) {
		Writer writer;
		write_common(writer, active_magic, value.common);
		writer.bytes(value.request_digest);
		writer.i64(value.created_at);
		writer.i64(value.last_seen);
		writer.u64(value.sequence);
		writer.vector(value.public_key);
		writer.string(value.key_reference);
		writer.string(value.fallback_reason);
		return writer.take();
	}

	Bytes serialize_active(const ActiveRecord& value) {
		Writer writer;
		const Bytes body = serialize_active_body(value);
		writer.bytes(body.data(), body.size());
		writer.vector(value.signature);
		return writer.take();
	}

	bool parse_active(const Bytes& data, ActiveRecord& value) {
		Reader reader(data);
		if (!read_common(reader, active_magic, value.common) || !reader.bytes(value.request_digest)
			|| !reader.i64(value.created_at) || !reader.i64(value.last_seen) || !reader.u64(value.sequence)
			|| !reader.vector(value.public_key) || !reader.string(value.key_reference)
			|| !reader.string(value.fallback_reason)) return false;
		value.signed_body.assign(data.begin(), data.begin() + static_cast<std::ptrdiff_t>(reader.offset()));
		return reader.vector(value.signature) && reader.finished();
	}

	bool parse_license(const Bytes& data, LicenseRecord& value) {
		Reader reader(data);
		std::array<std::uint8_t, 3> reserved{};
		if (!read_common(reader, license_magic, value.common) || !reader.bytes(value.request_digest)
			|| !reader.bytes(value.license_id) || !reader.i64(value.issued_at) || !reader.i64(value.expires_at)
			|| !reader.u64(value.features) || !reader.vector(value.public_key)
			|| !reader.u8(value.issuer_algorithm) || !reader.bytes(reserved.data(), reserved.size())) return false;
		if (reserved != std::array<std::uint8_t, 3>{} || value.issuer_algorithm != 1) return false;
		value.signed_body.assign(data.begin(), data.begin() + static_cast<std::ptrdiff_t>(reader.offset()));
		return reader.vector(value.signature, 1024) && reader.finished();
	}

	bool legacy_file(const Bytes& data) {
		return data.size() >= sizeof(legacy_magic) - 1
			&& std::equal(legacy_magic, legacy_magic + sizeof(legacy_magic) - 1, data.begin());
	}

	bool valid_binding_combination(const CommonRecord& record, const Bytes& public_key) {
		if (record.binding_type == pf::BindingType::Smbios)
			return record.key_algorithm == pf::DeviceKeyAlgorithm::None && public_key.empty();
		return record.binding_type == pf::BindingType::Tpm2Key
			&& record.key_algorithm != pf::DeviceKeyAlgorithm::None && !public_key.empty();
	}

	std::string normalize_field(std::string value, bool uuid) {
		const auto first = std::find_if_not(value.begin(), value.end(), [](unsigned char c) { return std::isspace(c) != 0; });
		const auto last = std::find_if_not(value.rbegin(), value.rend(), [](unsigned char c) { return std::isspace(c) != 0; }).base();
		value = first < last ? std::string(first, last) : std::string{};
		std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
		if (uuid) value.erase(std::remove_if(value.begin(), value.end(), [](char c) { return c == '-' || c == '{' || c == '}'; }), value.end());
		return value;
	}

	bool valid_smbios_field(const std::string& value, bool uuid) {
		if (value.empty()) return false;
		if (uuid && value.size() != 32) return false;
		std::string compact;
		for (unsigned char c : value) if (std::isalnum(c)) compact.push_back(static_cast<char>(c));
		if (!compact.empty() && (std::all_of(compact.begin(), compact.end(), [](char c) { return c == '0'; })
			|| std::all_of(compact.begin(), compact.end(), [](char c) { return c == 'F'; }))) return false;
		static const std::array<const char*, 7> invalid = { "UNKNOWN", "NONE", "DEFAULT STRING", "TO BE FILLED BY O.E.M.", "SYSTEM SERIAL NUMBER", "NOT SPECIFIED", "INVALID" };
		return std::none_of(invalid.begin(), invalid.end(), [&](const char* marker) { return value == marker; });
	}

#if defined(_WIN32)
	std::string smbios_string(const std::uint8_t* structure, std::size_t available, std::uint8_t index) {
		if (index == 0 || available < structure[1]) return {};
		const char* cursor = reinterpret_cast<const char*>(structure + structure[1]);
		const char* end = reinterpret_cast<const char*>(structure + available);
		for (std::uint8_t current = 1; cursor < end && *cursor != '\0'; ++current) {
			const char* terminator = std::find(cursor, end, '\0');
			if (current == index) return std::string(cursor, terminator);
			cursor = terminator == end ? end : terminator + 1;
		}
		return {};
	}

	bool read_smbios(std::string& uuid, std::string& board_serial, std::string& reason) {
		const DWORD signature = ('R' << 24) | ('S' << 16) | ('M' << 8) | 'B';
		const UINT size = GetSystemFirmwareTable(signature, 0, nullptr, 0);
		if (size < 8) { reason = "GetSystemFirmwareTable did not return SMBIOS data"; return false; }
		Bytes raw(size);
		if (GetSystemFirmwareTable(signature, 0, raw.data(), size) != size) { reason = "Unable to read SMBIOS data"; return false; }
		std::size_t offset = 8;
		while (offset + 4 <= raw.size()) {
			const auto* structure = raw.data() + offset;
			const std::uint8_t type = structure[0], length = structure[1];
			if (length < 4 || offset + length > raw.size()) break;
			std::size_t next = offset + length;
			while (next + 1 < raw.size() && !(raw[next] == 0 && raw[next + 1] == 0)) ++next;
			if (next + 1 >= raw.size()) break;
			const std::size_t available = next + 2 - offset;
			if (type == 1 && length >= 24) uuid = hex_string(structure + 8, 16);
			if (type == 2 && length >= 8) board_serial = smbios_string(structure, available, structure[7]);
			offset = next + 2;
			if (type == 127) break;
		}
		uuid = normalize_field(uuid, true);
		board_serial = normalize_field(board_serial, false);
		if (!valid_smbios_field(uuid, true) || !valid_smbios_field(board_serial, false)) {
			reason = "SMBIOS System UUID or baseboard serial is missing/invalid";
			return false;
		}
		return true;
	}
#else
	bool read_text_file(const char* path, std::string& value) {
		std::ifstream input(path);
		return input && static_cast<bool>(std::getline(input, value));
	}

	bool read_smbios(std::string& uuid, std::string& board_serial, std::string& reason) {
		if (!read_text_file("/sys/class/dmi/id/product_uuid", uuid)
			|| !read_text_file("/sys/class/dmi/id/board_serial", board_serial)) {
			reason = "Unable to read Linux DMI System UUID and baseboard serial";
			return false;
		}
		uuid = normalize_field(uuid, true);
		board_serial = normalize_field(board_serial, false);
		if (!valid_smbios_field(uuid, true) || !valid_smbios_field(board_serial, false)) {
			reason = "SMBIOS System UUID or baseboard serial is missing/invalid";
			return false;
		}
		return true;
	}
#endif

	bool create_smbios_binding(BindingMaterial& result, std::string& error) {
		std::string uuid, board_serial;
		if (!read_smbios(uuid, board_serial, error)) return false;
		Writer canonical;
		canonical.bytes(smbios_domain, sizeof(smbios_domain) - 1);
		canonical.string(uuid);
		canonical.string(board_serial);
		result.identity.binding_type = pf::BindingType::Smbios;
		result.identity.device_key_algorithm = pf::DeviceKeyAlgorithm::None;
		if (!sha256(canonical.data(), result.identity.machine_hash)) { error = "SHA-256 failed"; return false; }
		result.identity.display_code = display_code(result.identity.machine_hash);
		return true;
	}

#if defined(_WIN32)
	constexpr wchar_t tpm_ecdsa_key_name[] = L"MInDes.DeviceBinding.V2";
	constexpr wchar_t tpm_rsa_key_name[] = L"MInDes.DeviceBinding.V2.RSA";

	const wchar_t* windows_key_name(pf::DeviceKeyAlgorithm algorithm) {
		return algorithm == pf::DeviceKeyAlgorithm::Rsa2048 ? tpm_rsa_key_name : tpm_ecdsa_key_name;
	}

	std::string windows_key_reference(pf::DeviceKeyAlgorithm algorithm) {
		return algorithm == pf::DeviceKeyAlgorithm::Rsa2048 ? "MInDes.DeviceBinding.V2.RSA" : "MInDes.DeviceBinding.V2";
	}

	bool open_tpm_key(bool create, NCRYPT_PROV_HANDLE& provider, NCRYPT_KEY_HANDLE& key,
		pf::DeviceKeyAlgorithm& algorithm, std::string& error) {
		provider = 0; key = 0;
		if (NCryptOpenStorageProvider(&provider, MS_PLATFORM_CRYPTO_PROVIDER, 0) != ERROR_SUCCESS) {
			error = "Microsoft Platform Crypto Provider is unavailable"; return false;
		}
		for (const auto candidate : { pf::DeviceKeyAlgorithm::EcdsaP256, pf::DeviceKeyAlgorithm::Rsa2048 }) {
			if (NCryptOpenKey(provider, &key, windows_key_name(candidate), 0, 0) == ERROR_SUCCESS) {
				algorithm = candidate; return true;
			}
		}
		if (!create) { error = "MInDes TPM key is missing for the current user"; NCryptFreeObject(provider); provider = 0; return false; }

		SECURITY_STATUS status = NCryptCreatePersistedKey(provider, &key, NCRYPT_ECDSA_P256_ALGORITHM, tpm_ecdsa_key_name, 0, 0);
		if (status == ERROR_SUCCESS) status = NCryptFinalizeKey(key, 0);
		if (status == ERROR_SUCCESS) { algorithm = pf::DeviceKeyAlgorithm::EcdsaP256; return true; }
		if (key) { NCryptDeleteKey(key, 0); key = 0; }

		status = NCryptCreatePersistedKey(provider, &key, NCRYPT_RSA_ALGORITHM, tpm_rsa_key_name, 0, 0);
		DWORD rsa_bits = 2048;
		if (status == ERROR_SUCCESS) status = NCryptSetProperty(key, NCRYPT_LENGTH_PROPERTY,
			reinterpret_cast<PBYTE>(&rsa_bits), sizeof(rsa_bits), 0);
		if (status == ERROR_SUCCESS) status = NCryptFinalizeKey(key, 0);
		if (status == ERROR_SUCCESS) { algorithm = pf::DeviceKeyAlgorithm::Rsa2048; return true; }
		error = "Unable to create a current-user TPM ECDSA P-256 or RSA-2048 key";
		if (key) NCryptDeleteKey(key, 0);
		NCryptFreeObject(provider); provider = 0; key = 0; return false;
	}

	bool export_windows_tpm_public(NCRYPT_KEY_HANDLE key, pf::DeviceKeyAlgorithm algorithm,
		Bytes& public_key, std::string& error) {
		const wchar_t* blob_type = algorithm == pf::DeviceKeyAlgorithm::Rsa2048 ? BCRYPT_RSAPUBLIC_BLOB : BCRYPT_ECCPUBLIC_BLOB;
		DWORD size = 0;
		if (NCryptExportKey(key, 0, blob_type, nullptr, nullptr, 0, &size, 0) != ERROR_SUCCESS) { error = "Unable to size TPM public key"; return false; }
		Bytes blob(size);
		if (NCryptExportKey(key, 0, blob_type, nullptr, blob.data(), size, &size, 0) != ERROR_SUCCESS) { error = "Unable to export TPM public key"; return false; }
		if (algorithm == pf::DeviceKeyAlgorithm::EcdsaP256) {
			if (size < sizeof(BCRYPT_ECCKEY_BLOB) + 64) return false;
			const auto* header = reinterpret_cast<const BCRYPT_ECCKEY_BLOB*>(blob.data());
			if (header->dwMagic != BCRYPT_ECDSA_PUBLIC_P256_MAGIC || header->cbKey != 32) return false;
			public_key.resize(65); public_key[0] = 4;
			std::copy_n(blob.data() + sizeof(BCRYPT_ECCKEY_BLOB), 64, public_key.data() + 1);
			return true;
		}
		if (size < sizeof(BCRYPT_RSAKEY_BLOB)) return false;
		const auto* header = reinterpret_cast<const BCRYPT_RSAKEY_BLOB*>(blob.data());
		if (header->Magic != BCRYPT_RSAPUBLIC_MAGIC || header->BitLength != 2048 || header->cbModulus != 256 || header->cbPublicExp == 0) return false;
		Writer canonical;
		canonical.u32(header->cbPublicExp);
		canonical.bytes(blob.data() + sizeof(BCRYPT_RSAKEY_BLOB), header->cbPublicExp + header->cbModulus);
		public_key = canonical.take(); return true;
	}

	bool tpm_sign(const std::string& reference, const Hash256& digest, Bytes& signature, std::string& error) {
		NCRYPT_PROV_HANDLE provider = 0; NCRYPT_KEY_HANDLE key = 0;
		if (NCryptOpenStorageProvider(&provider, MS_PLATFORM_CRYPTO_PROVIDER, 0) != ERROR_SUCCESS) { error = "Microsoft Platform Crypto Provider is unavailable"; return false; }
		const bool rsa = reference == "MInDes.DeviceBinding.V2.RSA";
		if (NCryptOpenKey(provider, &key, rsa ? tpm_rsa_key_name : tpm_ecdsa_key_name, 0, 0) != ERROR_SUCCESS) { NCryptFreeObject(provider); error = "MInDes TPM key is missing"; return false; }
		BCRYPT_PKCS1_PADDING_INFO padding{ BCRYPT_SHA256_ALGORITHM };
		void* padding_info = rsa ? &padding : nullptr;
		DWORD flags = rsa ? NCRYPT_PAD_PKCS1_FLAG : 0, size = 0;
		SECURITY_STATUS status = NCryptSignHash(key, padding_info, const_cast<PBYTE>(digest.data()), static_cast<DWORD>(digest.size()), nullptr, 0, &size, flags);
		if (status == ERROR_SUCCESS) { signature.resize(size); status = NCryptSignHash(key, padding_info, const_cast<PBYTE>(digest.data()), static_cast<DWORD>(digest.size()), signature.data(), size, &size, flags); signature.resize(size); }
		NCryptFreeObject(key); NCryptFreeObject(provider);
		if (status != ERROR_SUCCESS) { error = "TPM signing failed"; signature.clear(); return false; }
		return true;
	}

	bool create_tpm_binding(bool create, BindingMaterial& result, std::string& error) {
#if defined(_DEBUG)
		if (GetEnvironmentVariableA("MINDES_TEST_DISABLE_TPM", nullptr, 0) != 0) {
			error = "TPM disabled by the Debug-only MINDES_TEST_DISABLE_TPM test hook";
			return false;
		}
#endif
		NCRYPT_PROV_HANDLE provider = 0; NCRYPT_KEY_HANDLE key = 0;
		pf::DeviceKeyAlgorithm algorithm = pf::DeviceKeyAlgorithm::None;
		if (!open_tpm_key(create, provider, key, algorithm, error)) return false;
		result.identity.binding_type = pf::BindingType::Tpm2Key;
		result.identity.device_key_algorithm = algorithm;
		const bool exported = export_windows_tpm_public(key, algorithm, result.identity.device_public_key, error);
		NCryptFreeObject(key); NCryptFreeObject(provider);
		if (!exported) return false;
		result.key_reference = windows_key_reference(algorithm);
		Writer canonical;
		canonical.bytes(tpm_domain, sizeof(tpm_domain) - 1);
		canonical.u8(static_cast<std::uint8_t>(result.identity.device_key_algorithm));
		canonical.vector(result.identity.device_public_key);
		if (!sha256(canonical.data(), result.identity.machine_hash)) { error = "SHA-256 failed"; return false; }
		result.identity.display_code = display_code(result.identity.machine_hash);
		return true;
	}
#else
	bool canonical_public_from_pem(const char* pem, pf::DeviceKeyAlgorithm& algorithm, Bytes& public_key) {
		BIO* bio = BIO_new_mem_buf(pem, -1);
		EVP_PKEY* pkey = bio ? PEM_read_bio_PUBKEY(bio, nullptr, nullptr, nullptr) : nullptr;
		bool ok = false;
		if (pkey && EVP_PKEY_base_id(pkey) == EVP_PKEY_EC) {
			EC_KEY* ec = EVP_PKEY_get1_EC_KEY(pkey);
			const EC_GROUP* group = ec ? EC_KEY_get0_group(ec) : nullptr;
			const EC_POINT* point = ec ? EC_KEY_get0_public_key(ec) : nullptr;
			public_key.resize(65);
			ok = group && point && EC_POINT_point2oct(group, point, POINT_CONVERSION_UNCOMPRESSED,
				public_key.data(), public_key.size(), nullptr) == 65;
			if (ec) EC_KEY_free(ec);
			if (ok) algorithm = pf::DeviceKeyAlgorithm::EcdsaP256;
		}
		else if (pkey && EVP_PKEY_base_id(pkey) == EVP_PKEY_RSA) {
			RSA* rsa = EVP_PKEY_get1_RSA(pkey); const BIGNUM *n = nullptr, *e = nullptr;
			if (rsa) RSA_get0_key(rsa, &n, &e, nullptr);
			const int modulus_size = n ? BN_num_bytes(n) : 0, exponent_size = e ? BN_num_bytes(e) : 0;
			if (modulus_size == 256 && exponent_size > 0 && exponent_size <= 8) {
				Bytes exponent(static_cast<std::size_t>(exponent_size)), modulus(256);
				if (BN_bn2binpad(e, exponent.data(), exponent_size) == exponent_size
					&& BN_bn2binpad(n, modulus.data(), 256) == 256) {
					Writer writer; writer.u32(static_cast<std::uint32_t>(exponent.size()));
					writer.bytes(exponent.data(), exponent.size()); writer.bytes(modulus.data(), modulus.size());
					public_key = writer.take(); algorithm = pf::DeviceKeyAlgorithm::Rsa2048; ok = true;
				}
			}
			if (rsa) RSA_free(rsa);
		}
		if (pkey) EVP_PKEY_free(pkey); if (bio) BIO_free(bio);
		if (!ok) public_key.clear(); return ok;
	}

	bool der_signature_to_raw(const std::uint8_t* der, std::size_t size, Bytes& raw) {
		const unsigned char* cursor = der;
		ECDSA_SIG* signature = d2i_ECDSA_SIG(nullptr, &cursor, static_cast<long>(size));
		if (!signature) return false;
		const BIGNUM *r = nullptr, *s = nullptr; ECDSA_SIG_get0(signature, &r, &s);
		raw.resize(64);
		const bool ok = BN_bn2binpad(r, raw.data(), 32) == 32 && BN_bn2binpad(s, raw.data() + 32, 32) == 32;
		ECDSA_SIG_free(signature); if (!ok) raw.clear(); return ok;
	}

#if defined(MINDES_HAS_TPM2_TSS)
	constexpr char linux_tpm_key_path[] = "/HS/SRK/mindes-device-v2";

	bool fapi_sign_digest(const Hash256& digest, Bytes& signature, Bytes* public_key,
		pf::DeviceKeyAlgorithm* key_algorithm, std::string& error, bool create) {
		FAPI_CONTEXT* context = nullptr;
		TSS2_RC status = Fapi_Initialize(&context, nullptr);
		if (status != TSS2_RC_SUCCESS) { error = "TPM2-TSS FAPI initialization failed"; return false; }
		if (create) Fapi_CreateKey(context, linux_tpm_key_path, "sign,noDa", "", "");
		uint8_t* der_signature = nullptr; size_t signature_size = 0; char* public_pem = nullptr; char* certificate = nullptr;
		status = Fapi_Sign(context, linux_tpm_key_path, nullptr, digest.data(), digest.size(), &der_signature, &signature_size, &public_pem, &certificate);
		Bytes canonical_public; pf::DeviceKeyAlgorithm algorithm = pf::DeviceKeyAlgorithm::None;
		bool ok = status == TSS2_RC_SUCCESS && public_pem
			&& canonical_public_from_pem(public_pem, algorithm, canonical_public);
		if (ok && algorithm == pf::DeviceKeyAlgorithm::EcdsaP256)
			ok = der_signature_to_raw(der_signature, signature_size, signature);
		else if (ok && algorithm == pf::DeviceKeyAlgorithm::Rsa2048)
			signature.assign(der_signature, der_signature + signature_size);
		if (ok && public_key) *public_key = std::move(canonical_public);
		if (ok && key_algorithm) *key_algorithm = algorithm;
		if (!ok) error = "TPM2-TSS could not create/open and use the current-user signing key";
		Fapi_Free(der_signature); Fapi_Free(public_pem); Fapi_Free(certificate); Fapi_Finalize(&context);
		return ok;
	}
#endif

	bool tpm_sign(const std::string&, const Hash256& digest, Bytes& signature, std::string& error) {
#if defined(MINDES_HAS_TPM2_TSS)
		return fapi_sign_digest(digest, signature, nullptr, nullptr, error, false);
#else
		error = "MInDes was built without TPM2-TSS support"; return false;
#endif
	}

	bool create_tpm_binding(bool create, BindingMaterial& result, std::string& error) {
#if defined(_DEBUG)
		if (std::getenv("MINDES_TEST_DISABLE_TPM") != nullptr) {
			error = "TPM disabled by the Debug-only MINDES_TEST_DISABLE_TPM test hook";
			return false;
		}
#endif
#if defined(MINDES_HAS_TPM2_TSS)
		Hash256 challenge{}; Bytes ignored;
		if (!random_bytes(challenge.data(), challenge.size()) || !fapi_sign_digest(challenge, ignored,
			&result.identity.device_public_key, &result.identity.device_key_algorithm, error, create)) return false;
		result.identity.binding_type = pf::BindingType::Tpm2Key;
		result.key_reference = linux_tpm_key_path;
		Writer canonical; canonical.bytes(tpm_domain, sizeof(tpm_domain) - 1); canonical.u8(static_cast<std::uint8_t>(result.identity.device_key_algorithm)); canonical.vector(result.identity.device_public_key);
		if (!sha256(canonical.data(), result.identity.machine_hash)) { error = "SHA-256 failed"; return false; }
		result.identity.display_code = display_code(result.identity.machine_hash); return true;
#else
		(void)create; (void)result; error = "MInDes was built without TPM2-TSS support"; return false;
#endif
	}
#endif

	bool get_enrollment_binding(BindingMaterial& result, std::string& error) {
		std::string tpm_error;
		if (create_tpm_binding(true, result, tpm_error)) return true;
		result = {};
		result.fallback_reason = tpm_error;
		if (!create_smbios_binding(result, error)) {
			error = "TPM unavailable (" + tpm_error + "); SMBIOS fallback failed (" + error + ")";
			return false;
		}
		result.fallback_reason = tpm_error;
		return true;
	}

	bool get_runtime_binding(pf::BindingType type, BindingMaterial& result, std::string& error) {
		if (type == pf::BindingType::Tpm2Key) return create_tpm_binding(false, result, error);
		if (type == pf::BindingType::Smbios) return create_smbios_binding(result, error);
		error = "Unknown binding type"; return false;
	}

	bool verify_ecdsa_p256(const Bytes& public_key, const Hash256& digest, const Bytes& signature) {
		if (public_key.size() != 65 || public_key[0] != 4 || signature.size() != 64) return false;
#if defined(_WIN32)
		BCRYPT_ALG_HANDLE algorithm = nullptr; BCRYPT_KEY_HANDLE key = nullptr;
		Bytes blob(sizeof(BCRYPT_ECCKEY_BLOB) + 64);
		auto* header = reinterpret_cast<BCRYPT_ECCKEY_BLOB*>(blob.data()); header->dwMagic = BCRYPT_ECDSA_PUBLIC_P256_MAGIC; header->cbKey = 32;
		std::copy_n(public_key.data() + 1, 64, blob.data() + sizeof(BCRYPT_ECCKEY_BLOB));
		bool ok = BCryptOpenAlgorithmProvider(&algorithm, BCRYPT_ECDSA_P256_ALGORITHM, nullptr, 0) == 0
			&& BCryptImportKeyPair(algorithm, nullptr, BCRYPT_ECCPUBLIC_BLOB, &key, blob.data(), static_cast<ULONG>(blob.size()), 0) == 0
			&& BCryptVerifySignature(key, nullptr, const_cast<PUCHAR>(digest.data()), static_cast<ULONG>(digest.size()), const_cast<PUCHAR>(signature.data()), static_cast<ULONG>(signature.size()), 0) == 0;
		if (key) BCryptDestroyKey(key); if (algorithm) BCryptCloseAlgorithmProvider(algorithm, 0); return ok;
#else
		EC_KEY* key = EC_KEY_new_by_curve_name(NID_X9_62_prime256v1); if (!key) return false;
		const unsigned char* cursor = public_key.data();
		if (!o2i_ECPublicKey(&key, &cursor, public_key.size())) { EC_KEY_free(key); return false; }
		ECDSA_SIG* sig = ECDSA_SIG_new(); BIGNUM* r = BN_bin2bn(signature.data(), 32, nullptr); BIGNUM* s = BN_bin2bn(signature.data() + 32, 32, nullptr);
		if (!sig || !r || !s || ECDSA_SIG_set0(sig, r, s) != 1) { if (sig) ECDSA_SIG_free(sig); else { BN_free(r); BN_free(s); } EC_KEY_free(key); return false; }
		const bool ok = ECDSA_do_verify(digest.data(), static_cast<int>(digest.size()), sig, key) == 1;
		ECDSA_SIG_free(sig); EC_KEY_free(key); return ok;
#endif
	}

	bool decode_rsa_public(const Bytes& public_key, Bytes& exponent, Bytes& modulus) {
		if (public_key.size() < 4 + 1 + 256) return false;
		const std::uint32_t exponent_size = static_cast<std::uint32_t>(public_key[0])
			| (static_cast<std::uint32_t>(public_key[1]) << 8)
			| (static_cast<std::uint32_t>(public_key[2]) << 16)
			| (static_cast<std::uint32_t>(public_key[3]) << 24);
		if (exponent_size == 0 || exponent_size > 8 || public_key.size() != 4 + exponent_size + 256) return false;
		exponent.assign(public_key.begin() + 4, public_key.begin() + 4 + exponent_size);
		modulus.assign(public_key.begin() + 4 + exponent_size, public_key.end());
		return true;
	}

	bool verify_rsa2048(const Bytes& public_key, const Hash256& digest, const Bytes& signature) {
		Bytes exponent, modulus;
		if (!decode_rsa_public(public_key, exponent, modulus) || signature.size() != 256) return false;
#if defined(_WIN32)
		BCRYPT_ALG_HANDLE algorithm = nullptr; BCRYPT_KEY_HANDLE key = nullptr;
		Bytes blob(sizeof(BCRYPT_RSAKEY_BLOB) + exponent.size() + modulus.size());
		auto* header = reinterpret_cast<BCRYPT_RSAKEY_BLOB*>(blob.data());
		header->Magic = BCRYPT_RSAPUBLIC_MAGIC; header->BitLength = 2048;
		header->cbPublicExp = static_cast<ULONG>(exponent.size()); header->cbModulus = 256;
		header->cbPrime1 = 0; header->cbPrime2 = 0;
		std::copy(exponent.begin(), exponent.end(), blob.begin() + sizeof(BCRYPT_RSAKEY_BLOB));
		std::copy(modulus.begin(), modulus.end(), blob.begin() + sizeof(BCRYPT_RSAKEY_BLOB) + exponent.size());
		BCRYPT_PKCS1_PADDING_INFO padding{ BCRYPT_SHA256_ALGORITHM };
		const bool ok = BCryptOpenAlgorithmProvider(&algorithm, BCRYPT_RSA_ALGORITHM, nullptr, 0) == 0
			&& BCryptImportKeyPair(algorithm, nullptr, BCRYPT_RSAPUBLIC_BLOB, &key, blob.data(), static_cast<ULONG>(blob.size()), 0) == 0
			&& BCryptVerifySignature(key, &padding, const_cast<PUCHAR>(digest.data()), static_cast<ULONG>(digest.size()),
				const_cast<PUCHAR>(signature.data()), static_cast<ULONG>(signature.size()), BCRYPT_PAD_PKCS1) == 0;
		if (key) BCryptDestroyKey(key); if (algorithm) BCryptCloseAlgorithmProvider(algorithm, 0); return ok;
#else
		RSA* rsa = RSA_new(); BIGNUM* n = BN_bin2bn(modulus.data(), static_cast<int>(modulus.size()), nullptr);
		BIGNUM* e = BN_bin2bn(exponent.data(), static_cast<int>(exponent.size()), nullptr);
		if (!rsa || !n || !e || RSA_set0_key(rsa, n, e, nullptr) != 1) {
			if (rsa) RSA_free(rsa); else { BN_free(n); BN_free(e); } return false;
		}
		const bool ok = RSA_verify(NID_sha256, digest.data(), static_cast<unsigned int>(digest.size()),
			signature.data(), static_cast<unsigned int>(signature.size()), rsa) == 1;
		RSA_free(rsa); return ok;
#endif
	}

	bool verify_device_digest(pf::DeviceKeyAlgorithm algorithm, const Bytes& public_key,
		const Hash256& digest, const Bytes& signature) {
		if (algorithm == pf::DeviceKeyAlgorithm::EcdsaP256) return verify_ecdsa_p256(public_key, digest, signature);
		if (algorithm == pf::DeviceKeyAlgorithm::Rsa2048) return verify_rsa2048(public_key, digest, signature);
		return false;
	}

	bool verify_device_signature(pf::DeviceKeyAlgorithm algorithm, const Bytes& public_key, const Bytes& body, const Bytes& signature) {
		Hash256 digest{}; if (!sha256(body, digest)) return false;
		return verify_device_digest(algorithm, public_key, digest, signature);
	}

	bool sign_device_body(const std::string& reference, const Bytes& body, Bytes& signature, std::string& error) {
		Hash256 digest{}; if (!sha256(body, digest)) { error = "SHA-256 failed"; return false; }
		return tpm_sign(reference, digest, signature, error);
	}

	bool verify_official_signature(const LicenseRecord& record, std::string& error) {
		if (!pf::license_key::configured) { error = "Official license public key is not configured"; return false; }
		Hash256 digest{}; if (!sha256(record.signed_body, digest)) { error = "SHA-256 failed"; return false; }
		Bytes public_key(65); public_key[0] = 4;
		std::copy(pf::license_key::ecdsa_p256_xy.begin(), pf::license_key::ecdsa_p256_xy.end(), public_key.begin() + 1);
		if (!verify_ecdsa_p256(public_key, digest, record.signature)) { error = "Official license signature is invalid"; return false; }
		return true;
	}

	bool equal_common(const CommonRecord& left, const CommonRecord& right) {
		return left.binding_type == right.binding_type && left.key_algorithm == right.key_algorithm
			&& left.machine_hash == right.machine_hash && left.request_id == right.request_id;
	}

	void report(bool enabled, const std::string& message, int color = 31) {
		if (!enabled) return;
		pf::printf_color_on_control("> (license) " + message, color);
		std::cout << std::endl;
	}
}

namespace pf {
	const char* binding_type_name(BindingType type) {
		switch (type) {
		case BindingType::Tpm2Key: return "TPM2_KEY";
		case BindingType::Smbios: return "SMBIOS";
		default: return "UNKNOWN";
		}
	}

	License& License::instance() { static License singleton; return singleton; }

	License::License()
		: executable_directory_(get_executable_directory()),
		license_path_(executable_directory_ / "mindes.license"),
		activate_path_(executable_directory_ / "mindes.activate"),
		request_path_(executable_directory_ / "mindes.request") {}

	std::filesystem::path License::resolve_path(const std::filesystem::path& supplied, const std::filesystem::path& default_path) const {
		if (supplied.empty()) return default_path;
		if (supplied.is_absolute()) return supplied.lexically_normal();
		return (executable_directory_ / supplied).lexically_normal();
	}

	ActivationResult License::activate_this_user(const std::filesystem::path& request_file) {
		ActivationResult result;
		result.activation_file = activate_path_;
		result.request_file = resolve_path(request_file, request_path_);
		is_mid_activated_ = false; is_license_init_ = false;

		BindingMaterial binding; std::string error;
		if (!get_enrollment_binding(binding, error)) { result.message = error; report(true, error); return result; }

		RequestRecord request;
		request.common.binding_type = binding.identity.binding_type;
		request.common.key_algorithm = binding.identity.device_key_algorithm;
		request.common.machine_hash = binding.identity.machine_hash;
		if (!random_bytes(request.common.request_id.data(), request.common.request_id.size())
			|| !random_bytes(request.nonce.data(), request.nonce.size())) { result.message = "Secure random generation failed"; report(true, result.message); return result; }
		request.created_at = unix_time_now(); request.public_key = binding.identity.device_public_key; request.fallback_reason = binding.fallback_reason;
		request.signed_body = serialize_request_body(request);
		if (request.common.binding_type == BindingType::Tpm2Key
			&& !sign_device_body(binding.key_reference, request.signed_body, request.proof, error)) { result.message = error; report(true, error); return result; }
		const Bytes request_data = serialize_request(request);

		ActiveRecord active;
		active.common = request.common; active.created_at = request.created_at; active.last_seen = request.created_at;
		active.public_key = request.public_key; active.key_reference = binding.key_reference; active.fallback_reason = binding.fallback_reason;
		if (!sha256(request_data, active.request_digest)) { result.message = "SHA-256 failed"; report(true, result.message); return result; }
		active.signed_body = serialize_active_body(active);
		if (active.common.binding_type == BindingType::Tpm2Key
			&& !sign_device_body(active.key_reference, active.signed_body, active.signature, error)) { result.message = error; report(true, error); return result; }

		if (!atomic_write(result.request_file, request_data, error) || !atomic_write(activate_path_, serialize_active(active), error)) {
			result.message = error; report(true, error); return result;
		}
		binding_type_ = binding.identity.binding_type; machine_code_ = binding.identity.display_code; is_mid_activated_ = true;
		result.success = true; result.binding_type = binding_type_; result.machine_code = machine_code_;
		result.message = binding.fallback_reason.empty() ? "Activation request generated with TPM 2.0 binding"
			: "TPM unavailable; activation request generated with SMBIOS binding: " + binding.fallback_reason;
		report(true, result.message, binding.fallback_reason.empty() ? 32 : 33);
		return result;
	}

	bool License::check_mid_active(bool debug) {
		is_mid_activated_ = false; is_license_init_ = false; binding_type_ = BindingType::Unknown; machine_code_.clear();
		Bytes data; if (!read_file(activate_path_, data)) { report(debug, "Unable to read activation file: " + activate_path_.string()); return false; }
		if (legacy_file(data)) { report(debug, "Legacy CPUID activation format is not supported; reactivate MInDes."); return false; }
		ActiveRecord active; if (!parse_active(data, active) || !valid_binding_combination(active.common, active.public_key)) { report(debug, "Invalid V2 activation file: " + activate_path_.string()); return false; }

		BindingMaterial current; std::string error;
		if (!get_runtime_binding(active.common.binding_type, current, error)) { report(debug, error); return false; }
		if (current.identity.device_key_algorithm != active.common.key_algorithm
			|| current.identity.machine_hash != active.common.machine_hash || current.identity.device_public_key != active.public_key) { report(debug, "Activation machine binding mismatch"); return false; }
		if (active.common.binding_type == BindingType::Tpm2Key) {
			if (!verify_device_signature(active.common.key_algorithm, active.public_key, active.signed_body, active.signature)) { report(debug, "Activation state signature is invalid"); return false; }
			Hash256 challenge{}; Bytes proof;
			if (!random_bytes(challenge.data(), challenge.size()) || !tpm_sign(active.key_reference, challenge, proof, error)
				|| !verify_device_digest(active.common.key_algorithm, active.public_key, challenge, proof)) { report(debug, "TPM private-key possession check failed: " + error); return false; }
		}
		const auto now = unix_time_now(); if (active.created_at > now || active.last_seen > now) { report(debug, "Local clock is earlier than activation state"); return false; }
		binding_type_ = active.common.binding_type; machine_code_ = display_code(active.common.machine_hash); is_mid_activated_ = true;
		return true;
	}

	bool License::check_license(bool debug, const std::filesystem::path& license_file) {
		is_license_init_ = false;
		if (!is_mid_activated_ && !check_mid_active(debug)) return false;
		Bytes active_data; ActiveRecord active;
		if (!read_file(activate_path_, active_data) || !parse_active(active_data, active)) return false;
		const auto path = resolve_path(license_file, license_path_);
		Bytes data; if (!read_file(path, data)) { report(debug, "Unable to read license file: " + path.string()); return false; }
		if (legacy_file(data)) { report(debug, "Legacy CPUID license format is not supported; reactivate MInDes."); return false; }
		LicenseRecord license; std::string error;
		if (!parse_license(data, license) || !valid_binding_combination(license.common, license.public_key)) { report(debug, "Invalid V2 license file: " + path.string()); return false; }
		if (!verify_official_signature(license, error)) { report(debug, error); return false; }
		if (!equal_common(active.common, license.common) || active.request_digest != license.request_digest || active.public_key != license.public_key) { report(debug, "License does not match this activation request"); return false; }
		const auto now = unix_time_now(); if (license.issued_at > now || license.expires_at <= now || license.expires_at <= license.issued_at) { report(debug, "License is not currently valid"); return false; }

		active.last_seen = now; ++active.sequence; active.signature.clear(); active.signed_body = serialize_active_body(active);
		if (active.common.binding_type == BindingType::Tpm2Key
			&& !sign_device_body(active.key_reference, active.signed_body, active.signature, error)) { report(debug, error); return false; }
		if (!atomic_write(activate_path_, serialize_active(active), error)) { report(debug, error); return false; }
		is_license_init_ = true; return true;
	}

	bool License::update_last_seen(bool debug) {
		if (!is_license_init_) {
			report(debug, "The license was not validated in this process; activation time was not updated");
			return false;
		}

		Bytes active_data;
		ActiveRecord active;
		if (!read_file(activate_path_, active_data)) {
			report(debug, "Unable to read activation file: " + activate_path_.string());
			return false;
		}
		if (!parse_active(active_data, active) || !valid_binding_combination(active.common, active.public_key)) {
			report(debug, "Invalid V2 activation file: " + activate_path_.string());
			return false;
		}

		BindingMaterial current;
		std::string error;
		if (!get_runtime_binding(active.common.binding_type, current, error)) {
			report(debug, error);
			return false;
		}
		if (current.identity.device_key_algorithm != active.common.key_algorithm
			|| current.identity.machine_hash != active.common.machine_hash
			|| current.identity.device_public_key != active.public_key) {
			report(debug, "Activation machine binding mismatch");
			return false;
		}
		if (active.common.binding_type == BindingType::Tpm2Key
			&& !verify_device_signature(active.common.key_algorithm, active.public_key, active.signed_body, active.signature)) {
			report(debug, "Activation state signature is invalid");
			return false;
		}

		const auto now = unix_time_now();
		if (active.created_at > now || active.last_seen > now) {
			report(debug, "Local clock is earlier than activation state");
			return false;
		}

		active.last_seen = now;
		++active.sequence;
		active.signature.clear();
		active.signed_body = serialize_active_body(active);
		if (active.common.binding_type == BindingType::Tpm2Key
			&& !sign_device_body(active.key_reference, active.signed_body, active.signature, error)) {
			report(debug, error);
			return false;
		}
		if (!atomic_write(activate_path_, serialize_active(active), error)) {
			report(debug, error);
			return false;
		}
		return true;
	}

	bool License::get_machine_identity(MachineIdentity& result) {
		Bytes data; ActiveRecord active;
		if (!read_file(activate_path_, data) || !parse_active(data, active)) return false;
		BindingMaterial binding; std::string error;
		if (!get_runtime_binding(active.common.binding_type, binding, error)
			|| binding.identity.machine_hash != active.common.machine_hash
			|| binding.identity.device_public_key != active.public_key) return false;
		result = binding.identity;
		return true;
	}

	std::string License::get_machine_code() { MachineIdentity identity; return get_machine_identity(identity) ? identity.display_code : std::string{}; }
	BindingType License::binding_type() const { return binding_type_; }
	bool License::is_license() const { return is_license_init_; }
	bool License::is_activated() const { return is_mid_activated_; }

	void License::print_activation_info() {
		Bytes data; ActiveRecord active;
		if (!read_file(activate_path_, data) || !parse_active(data, active)) { report(true, "MInDes is not activated"); return; }
		report(true, std::string("Binding type: ") + binding_type_name(active.common.binding_type), 36);
		report(true, "Machine code (128-bit display): " + display_code(active.common.machine_hash), 36);
		report(true, "Machine hash (256-bit): " + hex_string(active.common.machine_hash.data(), active.common.machine_hash.size()), 36);
		if (!active.fallback_reason.empty()) report(true, "TPM fallback reason: " + active.fallback_reason, 33);
	}

	void License::print_active_file_info() { print_activation_info(); }

	void License::print_license_info(const std::filesystem::path& license_file) {
		const auto path = resolve_path(license_file, license_path_); Bytes data; LicenseRecord license; std::string error;
		if (!read_file(path, data)) { report(true, "Unable to read license file: " + path.string()); return; }
		if (legacy_file(data)) { report(true, "Legacy CPUID license format is not supported; reactivate MInDes."); return; }
		if (!parse_license(data, license)) { report(true, "Invalid V2 license file: " + path.string()); return; }
		report(true, std::string("Binding type: ") + binding_type_name(license.common.binding_type), 36);
		report(true, "Machine code (128-bit display): " + display_code(license.common.machine_hash), 36);
		report(true, "Issued at (UTC epoch): " + std::to_string(license.issued_at), 36);
		report(true, "Expires at (UTC epoch): " + std::to_string(license.expires_at), 36);
		const bool signature_valid = verify_official_signature(license, error);
		report(true, signature_valid ? "Official signature: valid" : "Official signature: invalid (" + error + ")", signature_valid ? 32 : 31);
	}
}
