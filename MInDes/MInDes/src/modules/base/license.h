#pragma once

#include <array>
#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace pf {
	enum class BindingType : std::uint8_t {
		Unknown = 0,
		Tpm2Key = 1,
		Smbios = 2
	};

	enum class DeviceKeyAlgorithm : std::uint8_t {
		None = 0,
		EcdsaP256 = 1,
		Rsa2048 = 2
	};

	struct MachineIdentity {
		std::uint16_t binding_version = 2;
		BindingType binding_type = BindingType::Unknown;
		DeviceKeyAlgorithm device_key_algorithm = DeviceKeyAlgorithm::None;
		std::array<std::uint8_t, 32> machine_hash{};
		std::vector<std::uint8_t> device_public_key;
		std::string display_code;
	};

	struct ActivationResult {
		bool success = false;
		BindingType binding_type = BindingType::Unknown;
		std::filesystem::path activation_file;
		std::filesystem::path request_file;
		std::string machine_code;
		std::string message;
	};

	class License final {
	public:
		static License& instance();

		License(const License&) = delete;
		License& operator=(const License&) = delete;
		License(License&&) = delete;
		License& operator=(License&&) = delete;

		ActivationResult activate_this_user(
			const std::filesystem::path& request_file = {});
		bool check_mid_active(bool debug = true);
		bool check_license(
			bool debug = true,
			const std::filesystem::path& license_file = {});

		bool get_machine_identity(MachineIdentity& result);
		std::string get_machine_code();
		BindingType binding_type() const;
		bool is_license() const;
		bool is_activated() const;

		void print_activation_info();
		void print_active_file_info();
		void print_license_info(const std::filesystem::path& license_file = {});

	private:
		License();
		~License() = default;

		std::filesystem::path resolve_path(
			const std::filesystem::path& supplied,
			const std::filesystem::path& default_path) const;

		std::filesystem::path executable_directory_;
		std::filesystem::path license_path_;
		std::filesystem::path activate_path_;
		std::filesystem::path request_path_;
		BindingType binding_type_ = BindingType::Unknown;
		std::string machine_code_;
		bool is_mid_activated_ = false;
		bool is_license_init_ = false;
	};

	const char* binding_type_name(BindingType type);
}
